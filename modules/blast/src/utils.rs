use crate::consts::SCALE;
use genepred::GenePred;

/// Given a position within the concatenated exonic sequence, this function
/// maps it back to its corresponding genomic coordinate.
///
/// The function iterates through the gene's exons, calculating the length of
/// each exon. It keeps a running total of the length of the exons processed so
/// far (`current_pos`). If the target position (`pos`) falls within the current
/// exon, it calculates the offset from the exon's start and returns the absolute
/// genomic coordinate. If the target position is beyond all exons, it returns `None`.
///
/// # Arguments
///
/// * `pos` - A `u64` representing the position within the concatenated exonic sequence (0-indexed).
///
/// # Returns
///
/// An `Option<u64>` which is:
/// * `Some(genomic_position)` if the position is found within the exons.
/// * `None` if the position is outside the total length of all exons.
///
/// # Panics
///
/// This function does not panic under normal conditions.
///
fn get_pos_in_exons(record: &GenePred, pos: u64) -> Option<u64> {
    let mut current_pos = 0;
    let mut exons = record.exons();

    exons = match record.strand() {
        Some(genepred::Strand::Forward) => exons,
        Some(genepred::Strand::Reverse) => {
            // INFO: we scale the coordinates and then sort them
            // INFO: this is equivalent to reversing the exons and reversing the coordinates
            let mut exons = exons
                .iter()
                .map(|(start, end)| (SCALE - end, SCALE - start))
                .collect::<Vec<(u64, u64)>>();

            exons.sort_by(|a, b| a.0.cmp(&b.0));
            exons
        }
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    };

    for (exon_start, exon_end) in exons.iter() {
        let block_len = exon_end - exon_start; // INFO: exon length

        if pos < current_pos + block_len {
            // INFO: position falls inside this exon
            let offset = pos - current_pos;
            return Some(exon_start + offset);
        }

        current_pos += block_len;
    }

    // INFO: edge case -> dispute between 0-based and 1-based coords [always ends]
    if pos >= current_pos {
        return Some(exons.last().unwrap().1);
    }

    None
}

/// Maps the start and end positions of an Open Reading Frame (ORF) from the
/// concatenated exonic sequence to absolute genomic coordinates.
///
/// This function is crucial for converting ORF coordinates, which are relative to the
/// combined sequence of all exons, into the real genomic coordinates (chromosome, start,
/// end). It calls `get_pos_in_exons` for both the ORF's start and end positions.
///
/// The function also accounts for the gene's strand. For the forward strand, the
/// converted coordinates are returned directly. For the reverse strand, the coordinates
/// are flipped and adjusted based on a `SCALE` value to correctly reflect their position
/// on the reverse complement.
///
/// # Arguments
///
/// * `orf_start` - A `u64` representing the starting position of the ORF within the
///   concatenated exonic sequence.
/// * `orf_end` - A `u64` representing the ending position of the ORF within the
///   concatenated exonic sequence.
///
/// # Returns
///
/// A tuple `(u64, u64)` containing the absolute genomic start and end coordinates of the ORF.
///
/// # Panics
///
/// This function will panic if `get_pos_in_exons` returns `None` for either `orf_start`
/// or `orf_end`, indicating that the ORF coordinates are outside the exonic sequence.
///
pub fn get_cds_from_pos(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    let cds_start = get_pos_in_exons(record, orf_start).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_START to exonic coordinates: {orf_start} -> {record:?}")
    });

    let cds_end = get_pos_in_exons(record, orf_end).unwrap_or_else(|| {
        panic!("ERROR: could not map ORF_END to exonic coordinates: {orf_end} -> {record:?}")
    });

    match record.strand() {
        Some(genepred::Strand::Forward) => (cds_start, cds_end),
        Some(genepred::Strand::Reverse) => (SCALE - cds_end, SCALE - cds_start),
        Some(genepred::Strand::Unknown) | None => {
            panic!("ERROR: unexpected strand value: {:?}", record.strand())
        }
    }
}

/// Maps ORF transcript coordinates to absolute genomic coordinates (strand-aware).
///
/// This function takes ORF coordinates defined in **dataar transcript space**
/// (as produced by a CDS/ORF predictor) and returns the corresponding genomic
/// start and end coordinates of the coding region, accounting for exon structure
/// and strand orientation. For transcripts on the reverse strand, it assumes
/// coordinates have been **strand-normalized** using a fixed scale (`SCALE`) and
/// reverses them appropriately.
///
/// # Arguments
///
/// * `orf_start` - The start of the ORF in transcript (spliced) coordinates.
/// * `orf_end` - The end of the ORF in transcript (spliced) coordinates. This is exclusive.
///
/// # Returns
///
/// * `(u64, u64)` - A tuple of genomic start and end coordinates corresponding
///   to the CDS region. Returns `(start, end)` if mapping succeeds; otherwise,
///   may return `Error` if the coordinates do not map within the exon structure.
///
/// # Strand Behavior
///
/// * On the `Strand::Forward`, transcript coordinates map left-to-right across exons.
/// * On the `Strand::Reverse`, the exons are reversed and mapped using:
///   `SCALE - (exon_end - offset)`, where `SCALE` is a user-defined upper genomic bound
///   used for coordinate normalization.
///
/// # Example
///
/// ```rust, ignore
/// let mut gene = GenePred { /* initialized */ };
/// let orf_start = 15;
/// let orf_end = 45;
/// if let Some((cds_start, cds_end)) = gene.map_absolute_cds(orf_start, orf_end) {
///     println!("Genomic CDS coordinates: {} - {}", cds_start, cds_end);
/// }
/// ```
pub fn map_absolute_cds(record: &GenePred, orf_start: u64, orf_end: u64) -> (u64, u64) {
    // INFO: for Reverse this should be interpreted from 5' to 3'
    let (cds_start, cds_end) = get_cds_from_pos(record, orf_start, orf_end);

    assert!(
        cds_start > 0,
        "ERROR: {orf_start} is not bigger than 0! -> {record:?}"
    );
    assert!(
        cds_end > 0,
        "ERROR: {orf_end} is not bigget than 0! -> {record:?}"
    );

    assert!(
        cds_start >= record.start,
        "ERROR: {cds_start} is less than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.start
    );
    assert!(
        cds_end <= record.end,
        "ERROR: {cds_end} is greater than {}! Mapped ORF values {orf_start} - {orf_end} -> {record:?}",
        record.end
    );

    (cds_start, cds_end)
}
