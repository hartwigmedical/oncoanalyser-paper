WITH cohort AS (
  SELECT
    s._datasource,
    s.sampleId,
    s.donorId,
    MAX(NULLIF(TRIM(cs.primaryTumorLocation), '')) AS sample_tumor_location
  FROM sample AS s
  JOIN clinical_summary AS cs
    ON cs._datasource = s._datasource
   AND cs.sampleId = s.sampleId
  WHERE s.donorId IS NOT NULL
  GROUP BY s._datasource, s.sampleId, s.donorId
),

-- One tumour location per donor (most frequent across their samples, ties
-- broken alphabetically), so a donor is never counted in two location groups.
donor_location AS (
  SELECT
    _datasource,
    donorId,
    sample_tumor_location AS tumor_category
  FROM (
    SELECT
      _datasource,
      donorId,
      sample_tumor_location,
      ROW_NUMBER() OVER (
        PARTITION BY _datasource, donorId
        ORDER BY COUNT(*) DESC, sample_tumor_location
      ) AS location_rank
    FROM cohort
    WHERE sample_tumor_location IS NOT NULL
    GROUP BY _datasource, donorId, sample_tumor_location
  ) AS ranked
  WHERE location_rank = 1
),

t_cell_reads AS (
  SELECT
    _datasource,
    sampleId,
    SUM(assignedReads) AS t_cell_reads
  FROM cdr3Sequence
  WHERE locus IN ('TRA/TRD', 'TRB')
    AND filter = 'PASS'
  GROUP BY _datasource, sampleId
),

telomere AS (
  SELECT
    _datasource,
    sampleId,
    AVG(germlineTelomereLength) AS germline_telomere_length,
    AVG(somaticTelomereLength) AS somatic_telomere_length
  FROM telomereLength
  GROUP BY _datasource, sampleId
),

reported_virus AS (
  SELECT DISTINCT
    _datasource,
    sampleId,
    CASE
      WHEN UPPER(interpretation) = 'HPV'
        OR LOWER(virusName) LIKE '%papilloma%' THEN 'HPV'
      WHEN UPPER(interpretation) = 'EBV'
        OR LOWER(virusName) LIKE '%epstein%' THEN 'EBV'
      WHEN UPPER(interpretation) = 'MCV'
        OR LOWER(virusName) LIKE '%merkel%' THEN 'Merkel cell polyomavirus'
      WHEN UPPER(interpretation) = 'HBV'
        OR LOWER(virusName) LIKE '%hepatitis b%' THEN 'HBV'
      WHEN UPPER(interpretation) = 'HHV8'
        OR LOWER(virusName) LIKE '%herpesvirus 8%'
        OR LOWER(virusName) LIKE '%kaposi%' THEN 'HHV-8/KSHV'
      WHEN LOWER(virusName) LIKE '%cytomegalo%'
        OR LOWER(virusName) LIKE '%betaherpesvirus 5%' THEN 'CMV'
      ELSE 'Other virus'
    END AS virus_category
  FROM virusAnnotation
  WHERE reported = 1
),

virus_by_sample AS (
  SELECT
    _datasource,
    sampleId,
    CASE
      WHEN COUNT(DISTINCT virus_category) > 1 THEN 'Multiple viruses'
      ELSE MAX(virus_category)
    END AS virus_category
  FROM reported_virus
  GROUP BY _datasource, sampleId
),

rna_fusion_counts AS (
  SELECT
    _datasource,
    sampleId,
    COUNT(*) AS rna_fusion_count
  FROM rnaFusion
  GROUP BY _datasource, sampleId
)

SELECT
  c._datasource,
  c.sampleId,
  c.donorId,
  dl.tumor_category,
  CAST(COALESCE(tr.t_cell_reads, 0) AS DOUBLE PRECISION)
    / NULLIF(m.tumorMedianCoverage, 0) AS t_cell_score,
  t.germline_telomere_length,
  t.somatic_telomere_length,
  COALESCE(v.virus_category, 'No virus detected') AS virus_category,
  CASE
    WHEN rs.sampleId IS NULL THEN NULL
    WHEN UPPER(rs.qcStatus) LIKE '%FAIL%' THEN NULL
    ELSE COALESCE(rf.rna_fusion_count, 0)
  END AS rna_fusion_count
FROM cohort AS c
JOIN donor_location AS dl
  ON dl._datasource = c._datasource
 AND dl.donorId = c.donorId
LEFT JOIN t_cell_reads AS tr
  ON tr._datasource = c._datasource
 AND tr.sampleId = c.sampleId
LEFT JOIN metric AS m
  ON m._datasource = c._datasource
 AND m.sampleId = c.sampleId
LEFT JOIN telomere AS t
  ON t._datasource = c._datasource
 AND t.sampleId = c.sampleId
LEFT JOIN virus_by_sample AS v
  ON v._datasource = c._datasource
 AND v.sampleId = c.sampleId
LEFT JOIN rnaStatistics AS rs
  ON rs._datasource = c._datasource
 AND rs.sampleId = c.sampleId
LEFT JOIN rna_fusion_counts AS rf
  ON rf._datasource = c._datasource
 AND rf.sampleId = c.sampleId
ORDER BY c._datasource, c.donorId, c.sampleId
