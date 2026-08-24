import os
import shutil

from pie.annotation import (
    parse_annotations,
    _merge_intervals,
    GeneModel,
    _file_checksum,
    _load_or_create_db,
    _read_cached_checksum,
    _resolve_cache_path,
)


class TestParseAnnotations:
    def test_parse_gff3(self, gff3_file):
        genes = parse_annotations(gff3_file)
        assert len(genes) == 3

    def test_parse_gtf(self, gtf_file):
        genes = parse_annotations(gtf_file)
        assert len(genes) == 3

    def test_same_gene_count(self, gff3_file, gtf_file):
        assert len(parse_annotations(gff3_file)) == len(parse_annotations(gtf_file))


class TestGeneModel:
    def test_gene1_attributes(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g1 = genes[0]
        assert "gene1" in g1.gene_id.lower()
        assert g1.strand == "+"
        assert g1.cds_length == 90
        assert len(g1.cds_exons) == 1

    def test_gene2_multi_exon(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g2 = [g for g in genes if "gene2" in g.gene_id.lower()][0]
        assert len(g2.cds_exons) == 2
        assert g2.cds_length == 99  # 59 + 40
        assert g2.strand == "+"

    def test_gene3_minus_strand(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g3 = [g for g in genes if "gene3" in g.gene_id.lower()][0]
        assert g3.strand == "-"
        assert g3.cds_length == 81
        assert len(g3.cds_exons) == 1

    def test_cds_exons_sorted(self, gff3_file):
        genes = parse_annotations(gff3_file)
        for gene in genes:
            starts = [s for _, s, _ in gene.cds_exons]
            assert starts == sorted(starts)

    def test_coordinates_0based(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g1 = genes[0]
        # GFF says 1-90, should be 0-based: start=0, end=90 (half-open)
        assert g1.cds_exons[0] == ("chr1", 0, 90)


class TestReadOnlyFallback:
    """Annotation parsing must succeed when the GFF lives in a read-only directory."""

    def test_resolve_cache_path_returns_none_for_readonly_dir(self, gff3_file, tmp_path):
        """_resolve_cache_path returns None when the parent dir is not writable."""
        ro_dir = tmp_path / "readonly"
        ro_dir.mkdir()
        gff_copy = ro_dir / "genes.gff3"
        shutil.copy(gff3_file, gff_copy)
        os.chmod(ro_dir, 0o555)
        try:
            assert _resolve_cache_path(str(gff_copy)) is None
        finally:
            os.chmod(ro_dir, 0o755)

    def test_parse_annotations_readonly_dir(self, gff3_file, tmp_path):
        """parse_annotations succeeds (in-memory DB) when GFF dir is read-only."""
        ro_dir = tmp_path / "readonly"
        ro_dir.mkdir()
        gff_copy = ro_dir / "genes.gff3"
        shutil.copy(gff3_file, gff_copy)
        os.chmod(ro_dir, 0o555)
        try:
            genes = parse_annotations(str(gff_copy))
            assert len(genes) == 3
            # No .pie.db file should be created
            assert not (ro_dir / "genes.gff3.pie.db").exists()
        finally:
            os.chmod(ro_dir, 0o755)


class TestCacheReuse:
    """The cached .pie.db must actually be loaded on subsequent calls."""

    def test_cache_db_is_created(self, gff3_file, tmp_path):
        """First call creates the .pie.db file with a valid checksum."""
        gff_copy = tmp_path / "genes.gff3"
        shutil.copy(gff3_file, gff_copy)
        gff_str = str(gff_copy)

        parse_annotations(gff_str)

        cache_path = gff_str + ".pie.db"
        assert os.path.exists(cache_path)

        # The checksum stored in the DB must match the file's checksum
        import gffutils
        db = gffutils.FeatureDB(cache_path)
        stored = _read_cached_checksum(db)
        assert stored is not None
        assert stored == _file_checksum(gff_str)

    def test_cache_is_reused_on_second_call(self, gff3_file, tmp_path):
        """Second call with the same GFF reuses the cached DB (no rebuild)."""
        gff_copy = tmp_path / "genes.gff3"
        shutil.copy(gff3_file, gff_copy)
        gff_str = str(gff_copy)
        cache_path = gff_str + ".pie.db"

        # First call: creates cache
        genes1 = parse_annotations(gff_str)
        assert os.path.exists(cache_path)
        mtime1 = os.path.getmtime(cache_path)

        # Second call: should reuse — DB file mtime must not change
        genes2 = parse_annotations(gff_str)
        mtime2 = os.path.getmtime(cache_path)

        assert mtime1 == mtime2
        assert len(genes1) == len(genes2)
        for g1, g2 in zip(genes1, genes2):
            assert g1.gene_id == g2.gene_id
            assert g1.cds_exons == g2.cds_exons

    def test_cache_invalidated_when_gff_changes(self, gff3_file, gtf_file, tmp_path):
        """Cache is rebuilt when the GFF file is replaced with different content."""
        gff_copy = tmp_path / "genes.gff3"
        shutil.copy(gff3_file, gff_copy)
        gff_str = str(gff_copy)
        cache_path = gff_str + ".pie.db"

        # Build cache from original GFF
        parse_annotations(gff_str)
        checksum1 = _file_checksum(gff_str)

        # Overwrite with different content (GTF file — different checksum)
        shutil.copy(gtf_file, gff_copy)
        checksum2 = _file_checksum(gff_str)
        assert checksum1 != checksum2

        # Parse again — must rebuild, stored checksum must update
        genes = parse_annotations(gff_str)
        assert len(genes) == 3

        import gffutils
        db = gffutils.FeatureDB(cache_path)
        assert _read_cached_checksum(db) == checksum2


class TestCdsPhase:
    """The GFF phase column must be honoured on each transcript's first CDS.

    5'-partial models (no annotated start codon, contig-edge fragments) carry
    phase 1 or 2 there; ignoring it translates the whole transcript out of
    frame, so every codon's syn/nonsyn classification in that gene is wrong.
    """

    def _gff(self, tmp_path, name, lines):
        p = tmp_path / name
        p.write_text("##gff-version 3\n" + "\n".join(lines) + "\n")
        return str(p)

    def test_plus_strand_phase_trims_5_prime(self, tmp_path):
        gff = self._gff(tmp_path, "plus.gff", [
            "chr1\t.\tgene\t2\t91\t.\t+\t.\tID=g1",
            "chr1\t.\tmRNA\t2\t91\t.\t+\t.\tID=t1;Parent=g1",
            "chr1\t.\tCDS\t2\t91\t.\t+\t1\tID=c1;Parent=t1",
        ])
        (gene,) = parse_annotations(gff)
        # 0-based start 1 (GFF 2) + phase 1 -> 2; end unchanged
        assert gene.cds_exons == [("chr1", 2, 91)]
        assert gene.cds_length == 89

    def test_minus_strand_phase_trims_high_end(self, tmp_path):
        """On '-' the transcript's first CDS is the one with the highest end."""
        gff = self._gff(tmp_path, "minus.gff", [
            "chr1\t.\tgene\t2\t91\t.\t-\t.\tID=g2",
            "chr1\t.\tmRNA\t2\t91\t.\t-\t.\tID=t2;Parent=g2",
            "chr1\t.\tCDS\t2\t50\t.\t-\t0\tID=c2a;Parent=t2",
            "chr1\t.\tCDS\t60\t91\t.\t-\t2\tID=c2b;Parent=t2",
        ])
        (gene,) = parse_annotations(gff)
        assert gene.cds_exons == [("chr1", 1, 50), ("chr1", 59, 89)]

    def test_phase_zero_and_dot_leave_exons_alone(self, tmp_path):
        for i, phase in enumerate(("0", ".")):
            gff = self._gff(tmp_path, f"zero{i}.gff", [
                f"chr1\t.\tgene\t2\t91\t.\t+\t.\tID=g{i}",
                f"chr1\t.\tmRNA\t2\t91\t.\t+\t.\tID=t{i};Parent=g{i}",
                f"chr1\t.\tCDS\t2\t91\t.\t+\t{phase}\tID=c{i};Parent=t{i}",
            ])
            (gene,) = parse_annotations(gff)
            assert gene.cds_exons == [("chr1", 1, 91)], f"phase={phase!r}"

    def test_phase_spanning_a_short_first_exon(self, tmp_path):
        """A phase larger than the first exon carries over into the next one."""
        gff = self._gff(tmp_path, "short.gff", [
            "chr1\t.\tgene\t2\t91\t.\t+\t.\tID=g5",
            "chr1\t.\tmRNA\t2\t91\t.\t+\t.\tID=t5;Parent=g5",
            "chr1\t.\tCDS\t2\t2\t.\t+\t2\tID=c5a;Parent=t5",
            "chr1\t.\tCDS\t10\t91\t.\t+\t0\tID=c5b;Parent=t5",
        ])
        (gene,) = parse_annotations(gff)
        # 1 bp exon consumed entirely, 1 more base trimmed off the next
        assert gene.cds_exons == [("chr1", 10, 91)]


class TestFileChecksum:
    def test_edit_past_8kb_with_preserved_mtime_is_detected(self, tmp_path):
        """size+mtime+first-8KB let rsync -a / cp -p / touch -r hide an edit."""
        p = tmp_path / "big.gff"
        p.write_bytes(b"A" * 60000)
        before = _file_checksum(str(p))

        stat = os.stat(p)
        with open(p, "r+b") as f:
            f.seek(59611)
            f.write(b"B")  # same file size, well past the first 8 KB
        os.utime(p, ns=(stat.st_atime_ns, stat.st_mtime_ns))

        assert _file_checksum(str(p)) != before

    def test_mtime_change_alone_does_not_invalidate(self, tmp_path):
        """Content hash: touching a file must not force a DB rebuild."""
        p = tmp_path / "a.gff"
        p.write_bytes(b"A" * 100)
        before = _file_checksum(str(p))
        os.utime(p, ns=(0, 0))
        assert _file_checksum(str(p)) == before


class TestMergeIntervals:
    """Overlapping or repeated CDS features are merged when the exon list is
    built, so every consumer of cds_exons sees each base once."""

    def test_overlapping_merged(self):
        assert _merge_intervals([("c", 0, 10), ("c", 5, 20)]) == [("c", 0, 20)]

    def test_exact_duplicate_merged(self):
        assert _merge_intervals([("c", 0, 10), ("c", 0, 10)]) == [("c", 0, 10)]

    def test_disjoint_kept_separate(self):
        assert _merge_intervals([("c", 30, 40), ("c", 0, 10)]) == [
            ("c", 0, 10), ("c", 30, 40)]

    def test_different_contigs_never_merged(self):
        assert _merge_intervals([("c", 0, 10), ("d", 5, 20)]) == [
            ("c", 0, 10), ("d", 5, 20)]

    def test_gff_listing_a_cds_twice(self, tmp_path):
        """The whole point: a duplicated CDS must not double the CDS length."""
        p = tmp_path / "dup.gff"
        p.write_text("##gff-version 3\n" + "\n".join([
            "chr1\t.\tgene\t1\t90\t.\t+\t.\tID=g1",
            "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1",
            "chr1\t.\tCDS\t1\t90\t.\t+\t0\tID=c1;Parent=t1",
            "chr1\t.\tCDS\t1\t90\t.\t+\t0\tID=c2;Parent=t1",
        ]) + "\n")
        (gene,) = parse_annotations(str(p))
        assert gene.cds_exons == [("chr1", 0, 90)]
        assert gene.cds_length == 90

    def test_overlapping_cds_features(self, tmp_path):
        p = tmp_path / "ov.gff"
        p.write_text("##gff-version 3\n" + "\n".join([
            "chr1\t.\tgene\t1\t90\t.\t+\t.\tID=g1",
            "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1",
            "chr1\t.\tCDS\t1\t60\t.\t+\t0\tID=c1;Parent=t1",
            "chr1\t.\tCDS\t31\t90\t.\t+\t0\tID=c2;Parent=t1",
        ]) + "\n")
        (gene,) = parse_annotations(str(p))
        assert gene.cds_exons == [("chr1", 0, 90)]
