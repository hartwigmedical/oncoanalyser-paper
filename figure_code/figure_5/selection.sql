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
donor_location AS (
  SELECT _datasource, donorId, sample_tumor_location AS tumor_category
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
t_cell_by_sample AS (
  SELECT
    c._datasource,
    c.donorId,
    c.sampleId,
    CAST(COALESCE(SUM(cs.assignedReads), 0) AS DOUBLE PRECISION)
      / NULLIF(MAX(m.tumorMedianCoverage), 0) AS t_cell_score
  FROM cohort AS c
  LEFT JOIN cdr3Sequence AS cs
    ON cs._datasource = c._datasource
   AND cs.sampleId = c.sampleId
   AND cs.locus IN ('TRA/TRD', 'TRB')
   AND cs.filter = 'PASS'
  LEFT JOIN metric AS m
    ON m._datasource = c._datasource
   AND m.sampleId = c.sampleId
  GROUP BY c._datasource, c.donorId, c.sampleId
),
t_cell_by_donor AS (
  SELECT _datasource, donorId, AVG(t_cell_score) AS t_cell_score
  FROM t_cell_by_sample
  GROUP BY _datasource, donorId
),
telomere_by_donor AS (
  SELECT
    c._datasource,
    c.donorId,
    AVG(t.germlineTelomereLength) AS germline_telomere_length,
    AVG(t.somaticTelomereLength) AS somatic_telomere_length
  FROM cohort AS c
  LEFT JOIN telomereLength AS t
    ON t._datasource = c._datasource
   AND t.sampleId = c.sampleId
  GROUP BY c._datasource, c.donorId
),
reported_virus AS (
  SELECT DISTINCT
    c._datasource,
    c.donorId,
    CASE
      WHEN UPPER(COALESCE(v.interpretation, '')) = 'HPV'
        OR LOWER(COALESCE(v.virusName, '')) LIKE '%papilloma%' THEN 'HPV'
      WHEN UPPER(COALESCE(v.interpretation, '')) = 'EBV'
        OR LOWER(COALESCE(v.virusName, '')) LIKE '%epstein%' THEN 'EBV'
      WHEN UPPER(COALESCE(v.interpretation, '')) IN ('MCV', 'MCPYV')
        OR LOWER(COALESCE(v.virusName, '')) LIKE '%merkel%' THEN 'Merkel cell polyomavirus'
      WHEN UPPER(COALESCE(v.interpretation, '')) = 'HBV'
        OR LOWER(COALESCE(v.virusName, '')) LIKE '%hepatitis b%' THEN 'HBV'
      WHEN UPPER(COALESCE(v.interpretation, '')) IN ('HHV-8', 'KSHV') THEN 'HHV-8/KSHV'
      WHEN UPPER(COALESCE(v.interpretation, '')) = 'CMV'
        OR LOWER(COALESCE(v.virusName, '')) LIKE '%cytomegalo%' THEN 'CMV'
      ELSE 'Other virus'
    END AS virus_category
  FROM cohort AS c
  JOIN virusAnnotation AS v
    ON v._datasource = c._datasource
   AND v.sampleId = c.sampleId
  WHERE v.reported = 1
),
virus_by_donor AS (
  SELECT
    _datasource,
    donorId,
    CASE
      WHEN COUNT(DISTINCT virus_category) > 1 THEN 'Multiple viruses'
      ELSE MAX(virus_category)
    END AS virus_category
  FROM reported_virus
  GROUP BY _datasource, donorId
),
-- DNA fusion calls are mapped from their DNA sample to the donor before counting.
dna_fusions_by_donor AS (
  SELECT
    c._datasource,
    c.donorId,
    COUNT(DISTINCT sf.name) AS dna_fusion_count
  FROM cohort AS c
  JOIN svFusion AS sf
    ON sf._datasource = c._datasource
   AND sf.sampleId = c.sampleId
  GROUP BY c._datasource, c.donorId
),
usable_rna_donors AS (
  SELECT DISTINCT c._datasource, c.donorId
  FROM cohort AS c
  JOIN rnaStatistics AS rs
    ON rs._datasource = c._datasource
   AND rs.sampleId = c.sampleId
  WHERE UPPER(COALESCE(rs.qcStatus, '')) NOT LIKE '%FAIL_%'
),
rna_fusions_by_donor AS (
  SELECT
    c._datasource,
    c.donorId,
    COUNT(DISTINCT rf.name) AS rna_fusion_count
  FROM cohort AS c
  JOIN rnaStatistics AS rs
    ON rs._datasource = c._datasource
   AND rs.sampleId = c.sampleId
   AND UPPER(COALESCE(rs.qcStatus, '')) NOT LIKE '%FAIL_%'
  JOIN rnaFusion AS rf
    ON rf._datasource = c._datasource
   AND rf.sampleId = c.sampleId
  GROUP BY c._datasource, c.donorId
)
SELECT
  dl._datasource,
  dl.donorId,
  dl.tumor_category,
  tc.t_cell_score,
  tl.germline_telomere_length,
  tl.somatic_telomere_length,
  COALESCE(v.virus_category, 'No virus detected') AS virus_category,
  COALESCE(df.dna_fusion_count, 0) AS dna_fusion_count,
  CASE
    WHEN ur.donorId IS NULL THEN NULL
    ELSE COALESCE(rf.rna_fusion_count, 0)
  END AS rna_fusion_count
FROM donor_location AS dl
LEFT JOIN t_cell_by_donor AS tc
  ON tc._datasource = dl._datasource AND tc.donorId = dl.donorId
LEFT JOIN telomere_by_donor AS tl
  ON tl._datasource = dl._datasource AND tl.donorId = dl.donorId
LEFT JOIN virus_by_donor AS v
  ON v._datasource = dl._datasource AND v.donorId = dl.donorId
LEFT JOIN dna_fusions_by_donor AS df
  ON df._datasource = dl._datasource AND df.donorId = dl.donorId
LEFT JOIN usable_rna_donors AS ur
  ON ur._datasource = dl._datasource AND ur.donorId = dl.donorId
LEFT JOIN rna_fusions_by_donor AS rf
  ON rf._datasource = dl._datasource AND rf.donorId = dl.donorId
ORDER BY dl.donorId