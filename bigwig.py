"""Minimal pure-Python bigWig reader (region values only).

Written because pyBigWig has no win-64 wheel. Implements the bigWig spec
(Kent et al. 2010, Bioinformatics 26:2204): header, chromosome B+ tree,
R-tree index and the three section types (bedGraph, varStep, fixedStep).
Read-only; no zoom levels (full-resolution data are read for the query
interval, which is what the genome-track panels need).
"""
import struct, zlib
import numpy as np

BIGWIG_MAGIC = 0x888FFC26
CHROM_MAGIC = 0x78CA8C91
RTREE_MAGIC = 0x2468ACE0


class BigWig:
    def __init__(self, path):
        self.f = open(path, "rb")
        magic = struct.unpack("<I", self.f.read(4))[0]
        if magic == BIGWIG_MAGIC:
            self.e = "<"
        else:
            self.f.seek(0)
            magic = struct.unpack(">I", self.f.read(4))[0]
            if magic != BIGWIG_MAGIC:
                raise ValueError("not a bigWig file")
            self.e = ">"
        (self.version, self.zoomLevels, self.chromTreeOffset, self.fullDataOffset,
         self.fullIndexOffset, self.fieldCount, self.definedFieldCount,
         self.autoSqlOffset, self.totalSummaryOffset, self.uncompressBufSize,
         _res) = struct.unpack(self.e + "HHQQQHHQQIQ", self.f.read(60))
        self.chroms = self._read_chrom_tree()

    # ---------------------------------------------------------------- chrom B+
    def _read_chrom_tree(self):
        f, e = self.f, self.e
        f.seek(self.chromTreeOffset)
        magic, blockSize, keySize, valSize, itemCount, _ = struct.unpack(e + "IIIIQQ", f.read(32))
        if magic != CHROM_MAGIC:
            raise ValueError("bad chromosome B+ tree")
        chroms = {}

        def node(offset):
            f.seek(offset)
            isLeaf, _, count = struct.unpack(e + "BBH", f.read(4))
            if isLeaf:
                for _ in range(count):
                    key = f.read(keySize).rstrip(b"\x00").decode()
                    cid, csize = struct.unpack(e + "II", f.read(valSize))
                    chroms[key] = (cid, csize)
            else:
                kids = []
                for _ in range(count):
                    f.read(keySize)
                    kids.append(struct.unpack(e + "Q", f.read(8))[0])
                for k in kids:
                    node(k)

        node(self.chromTreeOffset + 32)
        return chroms

    # ------------------------------------------------------------------ R-tree
    def _overlapping_blocks(self, cid, start, end):
        f, e = self.f, self.e
        f.seek(self.fullIndexOffset)
        magic = struct.unpack(e + "I", f.read(4))[0]
        if magic != RTREE_MAGIC:
            raise ValueError("bad R-tree index")
        f.read(44)
        root = self.fullIndexOffset + 48
        out = []

        def cmp_lt(c1, b1, c2, b2):
            return (c1, b1) < (c2, b2)

        def node(offset):
            f.seek(offset)
            isLeaf, _, count = struct.unpack(e + "BBH", f.read(4))
            rows = f.read(count * (32 if isLeaf else 24))
            for i in range(count):
                if isLeaf:
                    sc, sb, ec, eb, off, size = struct.unpack_from(e + "IIIIQQ", rows, i * 32)
                    if not (cmp_lt(ec, eb, cid, start + 1) or cmp_lt(cid, end, sc, sb + 1)):
                        out.append((off, size))
                else:
                    sc, sb, ec, eb, off = struct.unpack_from(e + "IIIIQ", rows, i * 24)
                    if not (cmp_lt(ec, eb, cid, start + 1) or cmp_lt(cid, end, sc, sb + 1)):
                        node(off)

        node(root)
        return out

    # ------------------------------------------------------------------ values
    def intervals(self, chrom, start, end):
        """List of (start, end, value) overlapping [start, end)."""
        if chrom not in self.chroms:
            alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
            if alt not in self.chroms:
                return []
            chrom = alt
        cid = self.chroms[chrom][0]
        res = []
        for off, size in self._overlapping_blocks(cid, start, end):
            self.f.seek(off)
            buf = self.f.read(size)
            if self.uncompressBufSize > 0:
                buf = zlib.decompress(buf)
            bcid, bstart, bend, istep, ispan, btype, _, icount = struct.unpack_from(
                self.e + "IIIIIBBH", buf, 0)
            if bcid != cid:
                continue
            o = 24
            if btype == 1:      # bedGraph
                for _ in range(icount):
                    s, en, v = struct.unpack_from(self.e + "IIf", buf, o); o += 12
                    if en > start and s < end:
                        res.append((s, en, v))
            elif btype == 2:    # variable step
                for _ in range(icount):
                    s, v = struct.unpack_from(self.e + "If", buf, o); o += 8
                    en = s + ispan
                    if en > start and s < end:
                        res.append((s, en, v))
            elif btype == 3:    # fixed step
                s = bstart
                for _ in range(icount):
                    v = struct.unpack_from(self.e + "f", buf, o)[0]; o += 4
                    en = s + ispan
                    if en > start and s < end:
                        res.append((s, en, v))
                    s += istep
        res.sort()
        return res

    def values(self, chrom, start, end, nbins=None, agg="mean"):
        """Per-base array over [start, end), or binned means if nbins is given."""
        n = end - start
        a = np.full(n, np.nan, dtype=np.float32)
        for s, e_, v in self.intervals(chrom, start, end):
            a[max(0, s - start):min(n, e_ - start)] = v
        if nbins is None:
            return a
        idx = np.linspace(0, n, nbins + 1).astype(int)
        out = np.empty(nbins, dtype=np.float32)
        for i in range(nbins):
            seg = a[idx[i]:idx[i + 1]]
            seg = seg[~np.isnan(seg)]
            out[i] = (np.nanmean(seg) if agg == "mean" else np.nanmax(seg)) if seg.size else np.nan
        return out

    def close(self):
        self.f.close()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()
