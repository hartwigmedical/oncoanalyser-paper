-- GENETIC IMMUNE ESCAPE - LOH of HLA allele
WITH sample_hla AS (
  SELECT DISTINCT
    c._datasource,
    c.donorId,
    c.sampleId,
    NULLIF(TRIM(c.primaryTumorLocation), '') AS tumor_category,
    m.gieHlaLoh
  FROM clinical AS c
  JOIN markers AS m
    ON m._datasource = c._datasource
   AND m.sampleId = c.sampleId
  WHERE c.donorId IS NOT NULL
    AND m.gieHlaLoh IS NOT NULL
),
donor_location AS (
  SELECT
    _datasource,
    donorId,
    tumor_category
  FROM (
    SELECT
      _datasource,
      donorId,
      tumor_category,
      ROW_NUMBER() OVER (
        PARTITION BY _datasource, donorId
        ORDER BY COUNT(*) DESC, tumor_category
      ) AS location_rank
    FROM sample_hla
    WHERE tumor_category IS NOT NULL
    GROUP BY _datasource, donorId, tumor_category
  ) AS ranked
  WHERE location_rank = 1
),
donor_hla AS (
  SELECT
    _datasource,
    donorId,
    MAX(gieHlaLoh) AS has_hla_loh
  FROM sample_hla
  GROUP BY _datasource, donorId
)
SELECT
  d.donorId,
  l.tumor_category,
  d.has_hla_loh
FROM donor_hla AS d
JOIN donor_location AS l
  ON l._datasource = d._datasource
 AND l.donorId = d.donorId
ORDER BY l.tumor_category, d.donorId