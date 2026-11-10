WITH cohort AS (
  SELECT s._datasource,s.sampleId,s.donorId,MAX(NULLIF(TRIM(cs.primaryTumorLocation),'')) AS tumor_category
  FROM sample s JOIN clinical_summary cs ON cs._datasource=s._datasource AND cs.sampleId=s.sampleId
  WHERE s.donorId IS NOT NULL GROUP BY s._datasource,s.sampleId,s.donorId
), tc AS (
  SELECT _datasource,sampleId,SUM(assignedReads) AS reads FROM cdr3Sequence
  WHERE locus IN ('TRA/TRD','TRB') AND filter='PASS' GROUP BY _datasource,sampleId
), df AS (
  SELECT _datasource,sampleId,COUNT(*) AS dna_fusion_count FROM svFusion GROUP BY _datasource,sampleId
), rf AS (
  SELECT _datasource,sampleId,COUNT(*) AS rna_fusion_count FROM rnaFusion GROUP BY _datasource,sampleId
), rv AS (
  SELECT DISTINCT _datasource,sampleId,
    CASE WHEN UPPER(COALESCE(interpretation,''))='HPV' OR LOWER(COALESCE(virusName,'')) LIKE '%papilloma%' THEN 'HPV'
         WHEN UPPER(COALESCE(interpretation,''))='EBV' OR LOWER(COALESCE(virusName,'')) LIKE '%epstein%' THEN 'EBV'
         WHEN UPPER(COALESCE(interpretation,'')) IN ('MCV','MCPYV') OR LOWER(COALESCE(virusName,'')) LIKE '%merkel%' THEN 'Merkel cell polyomavirus'
         WHEN UPPER(COALESCE(interpretation,''))='HBV' OR LOWER(COALESCE(virusName,'')) LIKE '%hepatitis b%' THEN 'HBV'
         WHEN UPPER(COALESCE(interpretation,'')) IN ('HHV-8','KSHV') THEN 'HHV-8/KSHV'
         WHEN UPPER(COALESCE(interpretation,''))='CMV' OR LOWER(COALESCE(virusName,'')) LIKE '%cytomegalo%' THEN 'CMV'
         ELSE 'Other virus' END AS virus_category
  FROM virusAnnotation WHERE reported=1
), vb AS (
  SELECT _datasource,sampleId,CASE WHEN COUNT(DISTINCT virus_category)>1 THEN 'Multiple viruses' ELSE MAX(virus_category) END AS virus_category
  FROM rv GROUP BY _datasource,sampleId
)
SELECT c.sampleId,c.donorId,c.tumor_category,
  CAST(COALESCE(tc.reads,0) AS DOUBLE)/NULLIF(m.tumorMedianCoverage,0) AS t_cell_score,
  t.germlineTelomereLength AS germline_telomere_length,
  t.somaticTelomereLength AS somatic_telomere_length,
  COALESCE(vb.virus_category,'No virus detected') AS virus_category,
  COALESCE(df.dna_fusion_count,0) AS dna_fusion_count,
  CASE WHEN rs.sampleId IS NULL OR UPPER(COALESCE(rs.qcStatus,'')) LIKE '%FAIL_%' THEN NULL ELSE COALESCE(rf.rna_fusion_count,0) END AS rna_fusion_count
FROM cohort c
LEFT JOIN tc ON tc._datasource=c._datasource AND tc.sampleId=c.sampleId
LEFT JOIN metric m ON m._datasource=c._datasource AND m.sampleId=c.sampleId
LEFT JOIN telomereLength t ON t._datasource=c._datasource AND t.sampleId=c.sampleId
LEFT JOIN vb ON vb._datasource=c._datasource AND vb.sampleId=c.sampleId
LEFT JOIN df ON df._datasource=c._datasource AND df.sampleId=c.sampleId
LEFT JOIN rnaStatistics rs ON rs._datasource=c._datasource AND rs.sampleId=c.sampleId
LEFT JOIN rf ON rf._datasource=c._datasource AND rf.sampleId=c.sampleId
WHERE c.tumor_category IS NOT NULL
ORDER BY c.donorId,c.sampleId
