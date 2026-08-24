# Idea: M22 non-CRF claim-route after M21 CRF failure

- slug: `m22-noncrf-claim-route` · status: **untried** · added: 2026-06-23
- refs: DEC-001,M21-GENERANNO-1P2B-CRF-SCREEN,M19-GENERANNO-1P2B-RAWCAL-CLEANPLANTS

## Hypothesis
M21 refuted trained CRF decoding on GENERanno 1.2B, but M19 non-CRF still preserves a valuable low-FPR operating point. The next claim-bearing route should either improve non-CRF emissions/objective while preserving M19-like FPR, or transfer the FP-aware recipe to a cleaner-provenance backbone such as NT-v2.

## Why it matters
Prevents repeated CRF tuning after DEC-001 and makes the M22 local design gate explicit before any new GPU run.

## Next step
Compare two M22 branches locally: (A) non-CRF GENERanno objective/emission improvement as adaptation evidence; (B) clean-provenance NT-v2/backbone transfer for claim potential. Select one primary GPU direction and record stop criteria.

## Log
- 2026-06-23: created (status=untried)
