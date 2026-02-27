import os
import shutil

from pie.annotation import (
    parse_annotations,
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
