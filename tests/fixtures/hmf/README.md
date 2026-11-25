Synthetic HMF resource fixtures for tests/testthat/test-hmf.R, in the formats of the oncoanalyser
resource bundle (v2.x; test-hmf.R derives the v3 variants). Real gene names, **made-up** flags,
coordinates and transcripts: nothing here is HMF data. The first eight genes have four exons,
coordinates start..start+800 on their contig, and one non-canonical transcript that must be
ignored. TMPRSS2, ETV1, ERG and ETV4 exercise LINX's phasing: TMPRSS2 and ETV1 have 5'UTR
introns (ETV1 also a 3'UTR intron), ERG lies 0.9 Mb from TMPRSS2 on chr21.
