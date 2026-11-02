WITH cohort AS (
  SELECT s._datasource, s.sampleId, s.donorId,
         MAX(NULLIF(TRIM(cs.primaryTumorLocation), '')) AS tumor_category
  FROM sample s
  JOIN clinical_summary cs
    ON cs._datasource = s._datasource AND cs.sampleId = s.sampleId
  WHERE s.donorId IS NOT NULL
  GROUP BY s._datasource, s.sampleId, s.donorId
),
cdr3_counts AS (
  SELECT _datasource, sampleId,
         SUM(passSequences) AS cdr3_tcell_count
  FROM cdr3LocusSummary
  WHERE locus IN ('TRA/TRD', 'TRB')
  GROUP BY _datasource, sampleId
),
t_cell_reads AS (
  SELECT _datasource, sampleId,
         SUM(assignedReads) AS t_cell_reads
  FROM cdr3Sequence
  WHERE locus IN ('TRA/TRD', 'TRB') AND filter = 'PASS'
  GROUP BY _datasource, sampleId
),
telomere AS (
  SELECT _datasource, sampleId,
         AVG(germlineTelomereLength) AS germline_telomere_length,
         AVG(somaticTelomereLength) AS somatic_telomere_length
  FROM telomereLength
  GROUP BY _datasource, sampleId
),
telomere_rearrangement AS (
  SELECT DISTINCT _datasource, sampleId,
         1 AS has_telomere_rearrangement
  FROM structuralVariant
  WHERE LOWER(TRIM(COALESCE(insertSequenceRepeatClass, ''))) = 'satellite/telo'
),
reported_virus AS (
  SELECT DISTINCT _datasource, sampleId,
    CASE
      WHEN UPPER(COALESCE(interpretation, '')) = 'HPV'
        OR LOWER(COALESCE(virusName, '')) LIKE '%papilloma%' THEN 'HPV'
      WHEN UPPER(COALESCE(interpretation, '')) = 'EBV'
        OR LOWER(COALESCE(virusName, '')) LIKE '%epstein%' THEN 'EBV'
      WHEN UPPER(COALESCE(interpretation, '')) IN ('MCV', 'MCPYV')
        OR LOWER(COALESCE(virusName, '')) LIKE '%merkel%'
          THEN 'Merkel cell polyomavirus'
      WHEN UPPER(COALESCE(interpretation, '')) = 'HBV'
        OR LOWER(COALESCE(virusName, '')) LIKE '%hepatitis b%' THEN 'HBV'
      WHEN UPPER(COALESCE(interpretation, '')) IN ('HHV-8', 'KSHV')
          THEN 'HHV-8/KSHV'
      WHEN UPPER(COALESCE(interpretation, '')) = 'CMV'
        OR LOWER(COALESCE(virusName, '')) LIKE '%cytomegalo%'
        OR LOWER(COALESCE(virusName, '')) LIKE '%betaherpesvirus 5%' THEN 'CMV'
      ELSE 'Other virus'
    END AS virus_category
  FROM virusAnnotation
  WHERE reported = 1
),
virus_by_sample AS (
  SELECT _datasource, sampleId,
         CASE WHEN COUNT(DISTINCT virus_category) > 1
              THEN 'Multiple viruses'
              ELSE MAX(virus_category) END AS virus_category
  FROM reported_virus
  GROUP BY _datasource, sampleId
),
sample_metrics AS (
  SELECT c._datasource, c.sampleId, c.donorId,
         COALESCE(c.tumor_category, 'Unknown') AS tumor_category,
         COALESCE(cc.cdr3_tcell_count, 0) AS cdr3_tcell_count,
         COALESCE(cr.t_cell_reads, 0) AS t_cell_reads,
         m.tumorMedianCoverage AS sample_median_depth,
         COALESCE(cr.t_cell_reads, 0) /
           NULLIF(m.tumorMedianCoverage, 0) AS t_cell_score,
         t.germline_telomere_length,
         t.somatic_telomere_length,
         COALESCE(tr.has_telomere_rearrangement, 0)
           AS has_telomere_rearrangement,
         COALESCE(v.virus_category, 'No virus detected') AS virus_category
  FROM cohort c
  LEFT JOIN cdr3_counts cc
    ON cc._datasource = c._datasource AND cc.sampleId = c.sampleId
  LEFT JOIN t_cell_reads cr
    ON cr._datasource = c._datasource AND cr.sampleId = c.sampleId
  LEFT JOIN metric m
    ON m._datasource = c._datasource AND m.sampleId = c.sampleId
  LEFT JOIN telomere t
    ON t._datasource = c._datasource AND t.sampleId = c.sampleId
  LEFT JOIN telomere_rearrangement tr
    ON tr._datasource = c._datasource AND tr.sampleId = c.sampleId
  LEFT JOIN virus_by_sample v
    ON v._datasource = c._datasource AND v.sampleId = c.sampleId
),
patient_category AS (
  SELECT _datasource, donorId, tumor_category,
         AVG(cdr3_tcell_count) AS cdr3_tcell_count,
         AVG(t_cell_reads) AS t_cell_reads,
         AVG(sample_median_depth) AS sample_median_depth,
         AVG(t_cell_score) AS t_cell_score,
         AVG(germline_telomere_length) AS germline_telomere_length,
         AVG(somatic_telomere_length) AS somatic_telomere_length,
         MAX(has_telomere_rearrangement) AS has_telomere_rearrangement,
         CASE
           WHEN COUNT(DISTINCT CASE WHEN virus_category <> 'No virus detected'
                                    THEN virus_category END) = 0
             THEN 'No virus detected'
           WHEN COUNT(DISTINCT CASE WHEN virus_category <> 'No virus detected'
                                    THEN virus_category END) > 1
             THEN 'Multiple viruses'
           ELSE MAX(CASE WHEN virus_category <> 'No virus detected'
                         THEN virus_category END)
         END AS virus_category
  FROM sample_metrics
  GROUP BY _datasource, donorId, tumor_category
)
SELECT donorId, tumor_category, cdr3_tcell_count, t_cell_reads,
       sample_median_depth, t_cell_score AS "t-cell score",
       germline_telomere_length, somatic_telomere_length,
       has_telomere_rearrangement, virus_category
FROM patient_category
WHERE tumor_category <> 'Unknown'
ORDER BY donorId, tumor_category