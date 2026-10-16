# Complete CDS chain annotation Results working draft

2026-10-13 (Asia/Shanghai). Joint 24,576-base encoding exceeded the fixed Helixer comparator in complete CDS chain F1 on both species in the original TRAIN panel, but remained below ANNEVO and Tiberius. Combining that encoding with the native CDS prior did not preserve its gain. The completed three-objective DEV comparison also retained a tradeoff: structured continuation had the highest pooled chain F1, while the mixed objective performed better in rice but worse in Arabidopsis and pooled accuracy, with more reference-relative background chains. TRAIN fitting and historically reused DEV measurements do not establish independent generalization or publication readiness. The [Methods draft](M28_full_grammar_methods_working.md) defines the model and inference procedure.

## Initial complete chain accuracy on the TRAIN panel

We evaluated the final model after 4,608 updates on six TRAIN chromosomes from Arabidopsis thaliana and Oryza sativa, totaling 156,779,632 bases. Training used 1,536 fixed draws from 1,491 unique windows, repeated for three passes, rather than exhaustive supervision of all chains on those chromosomes. The three native tools received the same six FASTA sequences. Inference used no reference annotations, training labels, RNA or protein evidence. Evaluation selected the longest eligible CDS per assessable reference gene, with transcript-ID tie-breaking, and retained 17,739 Arabidopsis and 4,949 rice primary chains. References outside the model grammar remained in the denominator; native predictions were not filtered to the grammar. An exact match required the correct sequence identifier, strand and every ordered CDS interval. CDS phase consistency was assessed separately.

The initial model recovered 9,966 of the 22,688 primary reference chains among 21,032 unique predicted chains. Pooled precision was 47.385%, recall was 43.926% and F1 was 45.590%. ANNEVO, Tiberius and Helixer attained pooled F1 values of 75.542%, 72.655% and 63.592%, respectively, exceeding the initial model by 29.952, 27.065 and 18.002 percentage points. Pooled scores were calculated from pooled counts, rather than averaging chromosome scores.

| Method | Arabidopsis F1 % | Rice F1 % | Pooled F1 % | Species macro F1 % | Worst species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial grammar model | 50.487 | 28.278 | 45.590 | 39.383 | 28.278 |
| ANNEVO | 76.831 | 71.179 | 75.542 | 74.005 | 71.179 |
| Tiberius | 75.417 | 63.577 | 72.655 | 69.497 | 63.577 |
| Helixer | 68.122 | 50.294 | 63.592 | 59.208 | 50.294 |

The deficit was present in both species and was larger in rice. The initial model's rice F1 was 28.278%, compared with 50.487% in Arabidopsis; every comparator exceeded it in both species. Thus, the failure to match these tools was already evident within the original training range and cannot be attributed solely to transfer beyond that range. This comparison does not isolate a particular training objective, representation or structural module as its cause. The [completed comparison](../reports/M28-METHOD-RESTART/native_train_baseline_result.md) retains the underlying counts, per-species results and fixed operating points; the [benchmark figure](../reports/M28-METHOD-RESTART/figures/initial_train_benchmark/initial_train_complete_chain.svg) shows the same completed comparison.

## High CDS base accuracy did not imply accurate complete chains

The initial model achieved pooled CDS-base F1 of 90.360%, with base precision of 97.208% and recall of 84.414%, despite its substantially lower complete-chain F1. The corresponding CDS-base F1 values for ANNEVO, Tiberius and Helixer were 93.942%, 94.428% and 92.440%. This base metric measures the unstranded union of CDS intervals, whereas complete-chain matching requires all CDS boundaries and the strand to agree with one primary reference chain.

The initial model covered 0.618% of reference-relative background with predicted gene spans and produced 162 unique chains wholly within that background. The corresponding values were 1.190% and 383 for ANNEVO, 0.752% and 454 for Tiberius, and 3.613% and 1,614 for Helixer. Background was defined outside annotated gene-feature spans on either strand, not simply outside reference CDS intervals. Lower background burden and high base precision therefore did not establish complete-structure competitiveness.

All four methods had zero duplicate-chain records. Exact matches to other eligible reference isoforms, but not to the selected primary chain, numbered 748 for the initial model, 1,498 for ANNEVO, 1,552 for Tiberius and 1,339 for Helixer. These matches were reported separately rather than counted as primary true positives. A nonprimary match is not, by itself, evidence of a biologically false gene or a novel gene. These supplementary measures use the same [completed evaluation](../reports/M28-METHOD-RESTART/native_train_baseline_result.md).

## Completed continuation objectives on the TRAIN panel

We independently resumed the initial weights and saved AdamW state for a further 4,608 updates over three passes through the same ordered TRAIN draws. The backbone, labels, importance weights, optimizer settings and raw-FASTA caller remained fixed. The three branches retained the structured observation likelihood, replaced it with independent local role-observation likelihoods, or combined the structured and local negative log-likelihood terms with coefficients of 0.5 each. These coefficients do not equalize gradient magnitudes. We evaluated all three predeclared final weights at 9,216 updates, without intermediate-checkpoint selection, using the same six whole chromosomes and 22,688-reference evaluation denominator.

| Completed model | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Pooled CDS-base F1 % | Wholly background chains |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial structured model, 4,608 updates | 50.487 | 28.278 | 45.590 | 90.360 | 162 |
| Structured continuation, 9,216 updates | 61.836 | 35.413 | 55.745 | 93.200 | 197 |
| Local continuation, 9,216 updates | 44.160 | 33.382 | 41.729 | 94.097 | 491 |
| Mixed continuation, 9,216 updates | 51.673 | 38.327 | 48.628 | 94.158 | 712 |

Structured continuation correctly matched 12,696 primary chains among 22,862 unique predictions. Pooled precision was 55.533%, recall was 55.959% and F1 was 55.745%, exceeding the initial model by 10.155 percentage points. Arabidopsis and rice chain F1 increased by 11.349 and 7.135 percentage points, respectively. Species macro F1 was 48.624%, with rice remaining the worst species at 35.413%. All three continuations remained below the three fixed native comparators in both species; even the highest pooled continuation score was 7.847 percentage points below Helixer.

At the same update count, structured continuation exceeded local continuation in pooled chain F1 by 14.0164 percentage points, despite a 0.8965-percentage-point lower CDS-base F1. It produced 197 wholly background chains, compared with 491 for local continuation. Local continuation increased rice chain F1 relative to the initial model but reduced Arabidopsis and pooled chain F1, despite improving pooled CDS-base F1. These contrasts reinforce the distinction between CDS coverage and complete-chain recovery.

Mixed continuation correctly matched 11,212 primary chains among 23,425 unique predictions. Relative to structured continuation, its rice chain F1 was 2.9141 percentage points higher, but Arabidopsis and pooled F1 were 10.1625 and 7.1170 percentage points lower. Its species macro F1 was 45.000%, and its worst-species F1 was 38.327%. Mixed continuation had 0.9576-percentage-point higher pooled CDS-base F1 but 515 more wholly background chains; background-span fractions were 1.1435% for mixed and 0.2892% for structured continuation. The mixed objective therefore did not improve complete-chain accuracy uniformly across species or endpoints.

All three continuations retain the same complete-gene caller, so this is a contrast between training objectives rather than a decoder ablation. Matching update counts does not match gradient magnitudes or effective updates, and comparisons with the initial model also differ by 4,608 additional updates. These single-fit TRAIN diagnostics do not establish a structural module's necessity, independent generalization or a basis for changing the frozen DEV comparison. The [objective contrast record](../reports/M28-METHOD-RESTART/full_grammar_objective_contrast_contract.md#原-train-三目标完整比较) retains all three completed branches, their counts, background definitions and unchanged evaluation protocol.

## Native priors, sampling and encoding context on the TRAIN panel

Five further predeclared TRAIN branches completed their final 9,216-update weights, checkpoint and exposure verification, whole-chromosome raw inference and end-to-end evaluation. Each resumed the initial weights for 4,608 additional updates, without selecting an intermediate checkpoint. The native-prior branch added 28 zero-initialized coefficients over two frozen annotation-pretrained CDS heads. Fresh sampling used 4,608 fixed uniform draws from 4,219 unique windows; complete-core sampling used a fixed importance-corrected mixture with 4,107 unique windows. Joint encoding exposed each original 24,576-base read window to one native encoder forward rather than four independent 6,144-base blocks. The combined branch used joint encoding and both native CDS priors. All five retained the original six TRAIN chromosomes and 22,688 primary-reference denominator.

| Completed model or fixed comparator | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Wholly background chains |
| --- | ---: | ---: | ---: | ---: | ---: |
| Structured continuation | 61.836 | 35.413 | 55.745 | 48.624 | 197 |
| Native CDS prior | 69.198 | 47.980 | 64.339 | 58.589 | 259 |
| Fresh uniform draws | 67.288 | 48.869 | 63.073 | 58.078 | 498 |
| Joint 24,576-base encoding | 70.614 | 53.984 | 66.963 | 62.299 | 320 |
| Joint encoding plus native CDS prior | 68.009 | 43.042 | 62.451 | 55.525 | 190 |
| Complete-core mixture draws | 69.264 | 43.764 | 63.676 | 56.514 | 278 |
| ANNEVO | 76.831 | 71.179 | 75.542 | 74.005 | 383 |
| Tiberius | 75.417 | 63.577 | 72.655 | 69.497 | 454 |
| Helixer | 68.122 | 50.294 | 63.592 | 59.208 | 1,614 |

Joint encoding recovered 14,926 primary chains among 21,892 unique predictions, with pooled precision of 68.180%, recall of 65.788% and F1 of 66.963%. It exceeded structured continuation by 11.217 percentage points and the fixed Helixer comparator by 3.371 percentage points pooled, 2.492 points in Arabidopsis and 3.690 points in rice. It remained below ANNEVO and Tiberius in both species. Rice was the worst species for every method in this table, so the rice column also gives worst-species F1. These are TRAIN fitting comparisons, not evidence that the improvement transfers to an independent panel.

The combined prior-and-context branch attained pooled F1 of 62.451%, below joint encoding alone at 66.963% and the prior alone at 64.339%. Its 190 wholly background chains were fewer than either component branch's 320 and 259, respectively. Lower background burden therefore did not imply better complete-chain recovery. The prior-only and complete-core branches exceeded Helixer pooled but not in rice; fresh uniform sampling remained slightly below Helixer pooled. None exceeded ANNEVO or Tiberius in both species.

Matching update counts does not match effective optimization or upstream annotation exposure. The two native CDS heads explicitly reuse annotation-pretrained supervision; their inputs are not an original contribution of this model. Sampling changes finite exposure and gradient variation as well as the proposal, while the context contrast changes encoder visibility but keeps the whole-contig decoder. These single-fit descriptive comparisons do not establish a module's necessity or an additive causal mechanism. The [native-prior](../reports/M28-METHOD-RESTART/full_grammar_native_cds_prior_fit.md), [fresh-sampling](../reports/M28-METHOD-RESTART/full_grammar_fresh_draws_plan.md), [joint-encoding](../reports/M28-METHOD-RESTART/full_grammar_joint24_fit.md), [combined](../reports/M28-METHOD-RESTART/full_grammar_joint24_native_prior_plan.md) and [complete-core](../reports/M28-METHOD-RESTART/full_grammar_complete_core_draws_plan.md) records specify the fixed interventions. These results do not add the combined or complete-core models to the frozen DEV family.

## Completed sampling by encoding factorial comparison

The separately specified fresh-sampling joint-encoding branch completed its final-update-9,216 readout on the original six TRAIN chromosomes. All four caller FASTAs were directly equal to the normalized original DNA before the factorial comparison was interpreted; all four evaluations retained the same 22,688 primary references.

| Sampling and encoding | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % |
| --- | ---: | ---: | ---: |
| Original draws, split encoding | 61.836 | 35.413 | 55.745 |
| Fresh draws, split encoding | 67.288 | 48.869 | 63.073 |
| Original draws, joint encoding | 70.614 | 53.984 | 66.963 |
| Fresh draws, joint encoding | 73.368 | 59.636 | 70.303 |

Fresh-joint recovered 15,739 primary chains among 22,087 unique predictions, giving pooled precision of 71.259% and recall of 69.371%, species-macro F1 of 66.502% and worst-species F1 of 59.636%. Relative to original-draw joint encoding, primary hits increased by 813 net and predictions by 195; pooled F1 improved by 3.340 points, with gains of 2.754 points in Arabidopsis and 5.651 in rice. These aggregate differences do not establish individual recovered/lost-reference intersections.

The fixed F1-scale interaction F11−F10−F01+F00 was −3.988 percentage points pooled, −2.698 in Arabidopsis and −7.805 in rice. Both single modifications and their joint condition improved over the original split control, but the joint improvement was less than a simple sum on this nonlinear endpoint. A single-fit TRAIN interaction does not establish statistically significant antagonism, causal redundancy or module necessity; fresh sampling changes finite supervision exposure and gradient order, not coverage alone.

Fresh-joint exceeded the fixed Helixer operating point in both species but remained below Tiberius, ANNEVO and both default Carbon-A points. Its CDS-base F1 of 95.236% was slightly below original-joint's 95.390%, while reference-relative background-span coverage increased from 0.702% to 1.009% and wholly-background chains from 320 to 407. Better exact-chain F1 therefore did not mean every endpoint improved. The [completed four-cell result](../reports/M28-METHOD-RESTART/full_grammar_fresh_joint24_result.md) retains all chromosomes, species, operating points and negative contrasts. All B/C/D TRAIN accounting and their separate frozen six-method DEV comparison are now complete without candidate prescreening; the fixed-family results are reported below.

## Completed coarse information additional learning contrast

The informative and uniform-channel models completed the same original-TRAIN raw caller and common-ruler evaluation at final update 9,216. Both started after 4,608 informative updates; only the additional-fit epsilon changed from 0.1 to 2/3, retaining the fine-plus-coarse denominator and the constant coarse entropy term.

| Additional-fit channel | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Worst species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| Informative | 73.368 | 59.636 | 70.303 | 66.502 | 59.636 |
| Uniform | 74.276 | 58.099 | 70.645 | 66.187 | 58.099 |

Uniform recovered 15,846 primary chains among 22,173 predictions. Informative minus uniform was -162 hits and -0.908 F1 points in Arabidopsis, +55 hits and +1.537 points in rice, and -107 hits and -0.342 points pooled. All three chromosomes in each species had its corresponding direction. Pooled precision and recall differences were -0.206 and -0.472 points; macro F1 favored informative by 0.314 points despite lower pooled F1. These are aggregate net differences, not individual reference gain/loss intersections or independent replicate evidence.

Informative had slightly higher pooled CDS-base F1, 95.236% versus 95.111%, but higher background-span coverage, 1.009% versus 0.992%; wholly-background chains numbered 407 and 406. A verified additional learning signal therefore did not produce a consistent accuracy benefit across species and endpoints. The [completed contrast](../reports/M28-METHOD-RESTART/full_grammar_uniform_coarse_result.md) retains the unchanged 22,688-reference denominator and direct caller-DNA match. All three candidates remain in frozen B/C/D DEV regardless of TRAIN ranking. Their TRAIN accounting is complete; the initial DEV launch failed before data access because of a source-key naming mismatch. The schema-only R1 passed interface and real TRAIN-closure checks and completed all six raw tasks, three comparisons and the full six-method aggregation without changing the scientific protocol.

## Completed spliced CDS codon energy contrast on original TRAIN

Adding the fixed complete-spliced-CDS codon energy did not improve this matched additional-learning run. Both species and all six chromosomes had lower complete-chain F1 than fresh-joint at the final head-update-9,216 and codon-update-4,608 readout.

| Fixed TRAIN readout | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % |
| --- | ---: | ---: | ---: |
| Fresh-joint | 73.368 | 59.636 | 70.303 |
| Fresh-joint plus centered codon energy | 69.693 | 51.065 | 65.469 |

The codon branch recovered 14,976 primary chains among 23,062 predictions. It had 763 fewer primary hits and 975 more predictions than fresh-joint, with pooled F1 lower by 4.834 points, precision lower by 6.321 points and recall lower by 3.363 points. Species-macro F1 was 60.379%, and worst-species F1 was 51.065%. These aggregate net differences are not individual reference gain/loss intersections.

Pooled CDS-base F1 declined from 95.236% to 95.057%. Wholly-background chains increased from 407 to 551 even though background-span coverage declined from 1.009% to 0.990%; less covered background area did not imply fewer background structures or better complete-chain recovery. The [completed codon result](../reports/M28-METHOD-RESTART/full_grammar_codon_result.md) retains every chromosome, common denominator and negative contrast.

The verified zero-energy nesting, finite gradients and actual coding-carry implementation established integration, not accuracy. This single-fit intervention also changes parameter capacity, optimizer ages and jointly clipped trajectories; it neither isolates the codon term from the changed neural head nor rules out coding-content models generally. The negative branch remains in the frozen B/C/D DEV comparison without checkpoint selection, retraining or performance prescreening.


## Completed low rank encoder adaptation on original TRAIN

The fixed rank-four encoder-adaptation branch completed its final head-update-9,216 and adapter-update-4,608 readout, with the same starting head, fresh draws, partial observations, joint context and complete-set objective as fresh-joint. The intervention added 102,400 query/value residual parameters in the last four encoder layers; the original encoder weights stayed frozen, but its adapted representation did not.

| Fixed TRAIN readout | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % |
| --- | ---: | ---: | ---: |
| Fresh-joint, frozen encoder | 73.368 | 59.636 | 70.303 |
| Fresh-joint with rank-four q/v adaptation | 75.316 | 64.908 | 72.958 |
| Tiberius | 75.417 | 63.577 | 72.655 |
| ANNEVO | 76.831 | 71.179 | 75.542 |

Adaptation recovered 16,420 primary chains among 22,324 unique predictions, with pooled precision of 73.553% and recall of 72.373%. Relative to fresh-joint, pooled F1 increased by 2.656 points and primary hits by 681 net; F1 improved by 1.948 points in Arabidopsis and 5.272 in rice, with improvement on all six chromosomes. Species-macro F1 was 70.112%, and worst-species F1 was 64.908%. The net-hit difference alone does not identify individual recovery and loss intersections.

Adaptation exceeded Tiberius pooled F1 by 0.303 points and rice F1 by 1.331 points, but remained 0.101 points lower in Arabidopsis. It had 325 more pooled primary hits and 707 more predictions than Tiberius. It remained below ANNEVO by 2.583 points pooled and below both default Carbon-A operating points. Thus the competitive gap narrowed on this TRAIN panel without establishing uniform superiority or independently broader applicability.

Pooled CDS-base F1 increased from 95.236% to 95.435%, but CDS-base precision fell. Reference-relative background-span coverage rose from 1.009% to 1.090%, and wholly-background chains increased from 407 to 454. Improved complete-chain recovery therefore did not mean every output-quality endpoint improved. The [completed adaptation result](../reports/M28-METHOD-RESTART/full_grammar_joint_lora_result.md) preserves both species, all chromosomes, native comparators and background costs. Completed post-outcome reference replay separated 1,177 gained chains from 496 losses, with 15,243 retained hits and 5,772 missed by both models. Multisegment references contributed 639 of the 681 net additional hits (93.833%); they represented 76.221% of the reference population. Recall improved by 0.778 points for single-segment, 3.637 for two-to-five-segment and 3.773 for at-least-six-segment references. Among nine references spanning more than 24,576 bases, fresh-joint recovered one and adaptation two; only one net gain came from this small group. These fixed bins are descriptive recall accounting, not subset prediction precision/F1, independent long-gene capability or an isolated causal mechanism.

Added capacity, new adapter optimizer ages and altered jointly clipped optimization trajectories belong to the intervention. This single-fit TRAIN contrast does not isolate their individual causal effects, establish novelty or select a checkpoint from DEV. Both completed adaptation branches remain outside the frozen B/C/D DEV family.

## Completed full rank q v encoder adaptation on original TRAIN

The separate full-rank continuation completed its fixed final head-update-9,216 and native-q/v-update-4,608 readout. It resumed the same original head and AdamW moments as the frozen fresh-joint control, used the same 4,608 fresh draws from 4,219 unique windows, and retained all six original TRAIN contigs and 22,688 primary references. Only the eight query and value projection matrices in the last four encoder layers were adapted; all other encoder weights remained fixed. Direct caller-DNA equality was established before reference access.

| Method or fixed comparator | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Worst-species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| Frozen fresh-joint | 73.368 | 59.636 | 70.303 | 66.502 | 59.636 |
| Rank-four q/v | 75.316 | 64.908 | 72.958 | 70.112 | 64.908 |
| Full-rank q/v | 74.415 | 64.550 | 72.200 | 69.482 | 64.550 |
| ANNEVO | 76.831 | 71.179 | 75.542 | 74.005 | 71.179 |
| Tiberius | 75.417 | 63.577 | 72.655 | 69.497 | 63.577 |
| Helixer | 68.122 | 50.294 | 63.592 | 59.208 | 50.294 |

Full-rank adaptation recovered 16,106 primary chains among 21,927 predictions, giving pooled precision of 73.453%, recall of 70.989% and F1 of 72.200%. Relative to the frozen fresh-joint control, it added 367 net primary hits while making 160 fewer predictions. F1 increased on all six chromosomes, by 1.047 points in Arabidopsis, 4.914 in rice and 1.897 pooled. It remained 0.456 points below Tiberius pooled and 1.003 below it in Arabidopsis, although rice F1 was 0.974 higher. ANNEVO remained higher in both species and by 3.342 points pooled; the default native Carbon-A pipelines also remained higher on TRAIN.

Compared descriptively with the completed rank-four intervention, full-rank F1 was 0.901 points lower in Arabidopsis, 0.358 lower in rice and 0.758 lower pooled. It had 314 fewer primary hits and 397 fewer predictions pooled, but slightly higher CDS-base F1, 95.482% versus 95.435%, and fewer wholly-background chains, 359 versus 454. Background-span coverage was also lower, 0.886% versus 1.090%. This is an observed accuracy/background tradeoff between two fixed single-fit settings, not a rank-only causal experiment: parameterization, native learning rate, forward rounding and optimization trajectories also differ.

A post-outcome read-only replay of the unchanged TRAIN predictions retained 15,636 reference hits in both interventions. Full-rank adaptation gained 470 hits relative to rank-four and lost 784, with 5,798 missed by both. Thus its hit set was not a strict subset of the low-rank hit set. The fixed pre-existing segment bins accounted for 272 net lost hits in multisegment CDS chains and 42 in single-segment chains; the nine references spanning more than a native window had two low-rank hits and one full-rank hit. These intersections are exploratory reference-recall accounting, not subgroup precision/F1, a deployable ensemble, long-gene generalization or a validated rank mechanism.

The [completed full-rank result](../reports/M28-METHOD-RESTART/full_grammar_dense_qv_result.md) retains every chromosome, endpoint and negative comparator difference, along with source jobs and replay definitions. Both adaptation branches remain TRAIN-only and outside every frozen DEV family; no result selects a checkpoint or changes subsequent training, and sealed data remain closed.

## Native learning rate had opposite effects in low rank and full rank adaptation

The completed fixed four-cell TRAIN contrast filled the two missing native learning-rate settings while keeping the original rank-four high-rate and full-rank low-rate outputs unchanged. Every cell used one fixed fit and seed, the same starting update-4,608 head, ordered draws and observations, and the sole predefined head-update-9,216/native-update-4,608 readout. All 22,688 primary references remained in the denominator.

| Native q/v setting | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Worst species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rank four, native lr 10⁻⁵ | 70.415 | 51.449 | 65.929 | 60.932 | 51.449 |
| Rank four, native lr 3 × 10⁻⁴ | 75.316 | 64.908 | 72.958 | 70.112 | 64.908 |
| Full rank, native lr 10⁻⁵ | 74.415 | 64.550 | 72.200 | 69.482 | 64.550 |
| Full rank, native lr 3 × 10⁻⁴ | 65.132 | 37.972 | 58.982 | 51.552 | 37.972 |

Increasing the native rate improved rank-four pooled chain F1 by 7.030 percentage points but reduced full-rank F1 by 13.218 points. The corresponding Arabidopsis/rice changes were +4.901/+13.458 and −9.283/−26.578 points, respectively. The predefined difference in differences, full-rank high-minus-low minus rank-four high-minus-low, was −20.247 points pooled. At the lower matched rate, full rank exceeded rank four by 6.271 pooled points; at the higher rate, rank four exceeded full rank by 13.976. These opposite sensitivities preclude reading the original differently tuned E/M comparison as a universal rank advantage.

The new rank-four low-rate cell recovered 15,720 primary chains among 25,000 predictions, with 1,092 wholly-background chains. The new full-rank high-rate cell recovered 13,412 among 22,790, with 428 wholly-background chains and CDS-base F1 of 95.361% despite chain F1 of 58.982%. Neither new cell improved its parameterization's original tested setting. Native-rate matching still does not match parameterization, forward rounding, joint clipping or optimizer trajectories, and each cell has only one fit/seed. The [complete four-cell result](../reports/M28-METHOD-RESTART/full_grammar_qv_lr_controls_result.md) retains all chromosome, species, background and negative contrasts. This TRAIN-only optimization-sensitivity readout selects no independent-validation winner, adds no frozen DEV method and does not establish broader applicability or competitive superiority.

## Strict E M agreement on completed TRAIN calls

After the fixed rank-four and full-rank outcomes were known, a separate TRAIN-only contrast retained exactly the CDS geometries predicted by both models. Agreement required the same sequence identifier, strand and every ordered CDS interval. The rule was fixed before this readout and did not use reference matches, role energies or posterior thresholds. It is a post-outcome diagnostic of completed raw calls, not a preregistered independent validation or an implemented new raw-FASTA caller. Original weights and predictions were unchanged, and the complete parent metrics and paired DNA inputs were replayed before interpretation.

The intersection retained 19,221 unique chains and 15,636 primary matches among all 22,688 references. Pooled precision, recall and F1 were 81.349%, 68.917% and 74.619%, respectively. F1 exceeded rank-four adaptation by 1.660 percentage points and full-rank adaptation by 2.419 points; both species improved relative to both parents. It also exceeded the fixed Tiberius endpoint in each species and pooled accuracy, but remained below ANNEVO.

| Completed TRAIN readout | Arabidopsis F1 % | Rice F1 % | Pooled F1 % |
| --- | ---: | ---: | ---: |
| Rank-four q/v | 75.316 | 64.908 | 72.958 |
| Full-rank q/v | 74.415 | 64.550 | 72.200 |
| Exact E/M intersection | 76.241 | 68.540 | 74.619 |
| Tiberius | 75.417 | 63.577 | 72.655 |
| ANNEVO | 76.831 | 71.179 | 75.542 |

The higher complete-chain F1 did not represent a uniform accuracy gain. Relative to the rank-four output, the rule removed 3,103 predictions, including 784 primary matches, reducing complete-chain recall by 3.456 points. Relative to full-rank output, it removed 2,706 predictions and 470 primary matches, reducing recall by 2.072 points. Pooled CDS-base F1 fell to 90.389%, from 95.435% and 95.482% for the two parents. Wholly-background predictions declined to 245 and background span burden to 0.401%, but these reference-relative reductions do not identify every discarded prediction as biologically false.

The remaining pooled F1 gap to ANNEVO was 0.923 points, with gaps of 0.590 in Arabidopsis and 2.639 in rice. Complete-chain agreement thus supplied an exploratory precision–recall tradeoff, not independently established superiority or broader applicability. The [fixed-rule result and full accounting](../reports/M28-METHOD-RESTART/full_grammar_e_m_consensus_train.md) preserve all chromosomes, source checks, counts, isoform categories and costs. The selection-only CPU time excludes both parents' fitting and raw inference, and no frozen DEV family was expanded.


### Equal E M energy pooling on completed TRAIN calls

A separate fixed contrast averaged the two models' aligned fifteen-role logits with weights 0.5 and 0.5, promoting the original FP32 scores to FP64 before addition. It then solved the unchanged complete gene-set MAP problem over each entire oriented TRAIN chromosome. For a shared grammar this averages complete-set energies and is equivalent to an equal logarithmic opinion pool; it is not an intersection, union selector or average of gene probabilities. No weight, temperature or threshold was selected. Both parent outcomes and the hard-agreement result were already known, making this a post-outcome TRAIN exploration rather than independent validation.

The energy pool produced 22,066 unique chains with 16,527 primary matches among the unchanged 22,688 references. Pooled precision, recall and F1 were 74.898%, 72.845% and 73.857%, respectively. Arabidopsis F1 was 75.439% and rice F1 68.413%; species macro and worst-species F1 were 71.926% and 68.413%. Both species were numerically above the two parents and Tiberius, but the Arabidopsis difference from Tiberius was only 0.022 points. Pooled F1 remained 1.685 points below ANNEVO and 0.762 below hard agreement.

| Fixed TRAIN readout | Pooled precision % | Recall % | Chain F1 % | CDS-base F1 % |
| --- | ---: | ---: | ---: | ---: |
| Rank-four q/v | 73.553 | 72.373 | 72.958 | 95.435 |
| Full-rank q/v | 73.453 | 70.989 | 72.200 | 95.482 |
| Exact E/M intersection | 81.349 | 68.917 | 74.619 | 90.389 |
| Equal E/M energy pool | 74.898 | 72.845 | 73.857 | 95.513 |
| ANNEVO | 78.740 | 72.593 | 75.542 | 93.942 |

Compared with hard agreement, the energy pool recovered 891 additional primary matches and raised recall by 3.927 points and CDS-base F1 by 5.124 points. Precision fell by 6.450 points, yielding lower complete-chain F1. Its 379 wholly-background predictions and 0.924% background-span coverage were lower than rank-four output but higher than full-rank output; background span also exceeded Tiberius and hard agreement. Restoring coverage therefore did not remove the precision or background tradeoff.

The pool retained all 15,636 hard-agreement primary hits. Its 891 additional hits consisted of 518 CDS geometries previously emitted only by rank-four adaptation, 342 emitted only by full-rank adaptation and 31 emitted by neither. Of 462 geometries absent from both parents' outputs, 31 matched a primary reference and five matched another assessable isoform. The remaining 426 were neither-exact benchmark residuals, not proven false genes. “New” here refers to output geometry rather than biological gene discovery, and most additional primary hits came from retaining one parent's existing chain.

The [completed fixed-rule result](../reports/M28-METHOD-RESTART/full_grammar_e_m_energy_pool_train.md) retains all chromosomes, reference strata, endpoint and background contrasts, geometry partitions and CPU costs. Its whole-TRAIN readout used already-completed raw role arrays, so CPU redecoding time excludes parent fitting and native encoding. A separate [strict raw-FASTA caller qualification](../reports/M28-METHOD-RESTART/full_grammar_e_m_raw_native.md) independently loaded both final checkpoints and reproduced the equal-energy decoder and GFF outputs on two fixed 49,152-base TRAIN views and a short artificial sequence. Measured peak allocated VRAM was 4.536 GiB. This established input-only engineering integration, not fresh native reencoding and accuracy evaluation of all six chromosomes. No frozen DEV family was expanded, and these TRAIN endpoints do not establish broader applicability, statistical superiority or publication readiness.

A separate [reference geometry replay](../reports/M28-METHOD-RESTART/full_grammar_e_m_geometry_train.md) retained every primary reference and reproduced the unchanged four prediction sets. Among the energy pool's 6,161 missed primary chains, 1,855 (30.109%) had no overlapping same-strand prediction span, 1,632 (26.489%) had both correct outer CDS boundaries in one prediction but a different internal structure, and 1,482 (24.055%) had the correct transcriptional stop boundary but not the start boundary. These descriptive categories do not isolate a causal start-codon or splice mechanism. Pooling gained 373 and lost 266 primary hits relative to rank-four adaptation; its Arabidopsis TP count fell by 31 while rice TP rose by 138. The Arabidopsis F1 increase therefore cannot be described as a recall increase.


## Fixed local role disagreement penalty reduced TRAIN chain F1

A separate fixed TRAIN contrast replaced the mean centered role energy with the lower of the two branches' centered energies at every base and role. Equivalently, it subtracted one half of their absolute disagreement from the mean. Both parent outcomes were already known; no penalty coefficient, threshold or temperature was selected. The grammar, complete-set decoder, original weights and full reference denominator remained unchanged.

The rule recovered 15,990 primary chains among 21,257 predictions. Pooled precision, recall and F1 were 75.222%, 70.478% and 72.773%; species macro and worst-species F1 were 70.371% and 66.040%. Arabidopsis F1 was 74.701% and rice F1 66.040%, respectively 0.738 and 2.373 percentage points below mean energy pooling. Relative to that pool, it emitted 809 fewer chains but recovered 537 fewer primary references: a 0.324-point precision increase accompanied a 2.367-point recall loss and a 1.084-point F1 decrease. Pooled CDS-base F1 was 95.482%, slightly below the mean pool's 95.513%.

Background-span coverage fell from 0.924% to 0.671%, and wholly-background chains from 379 to 264, but these improvements did not compensate for complete-chain losses. The new output retained 15,367 of the 15,636 exact hits shared by both parents. It also emitted 759 geometries absent from both parents, of which 92 matched a primary reference, 82 another assessable isoform and 585 neither; these categories do not imply novel or biologically false genes. Pooled F1 remained 2.769 points below ANNEVO and 5.399 below the fixed filtered Carbon-A comparator, with deficits in both species. This fixed penalty therefore did not improve the mean pool's complete-chain accuracy on the original TRAIN panel. The [complete negative result](../reports/M28-METHOD-RESTART/full_grammar_minimum_centered_energy_train_plan.md) preserves every chromosome, comparator and reference stratum. Its cached-role redecoding is separate from the [finite raw-FASTA native qualification](../reports/M28-METHOD-RESTART/full_grammar_minimum_centered_raw.md), and neither constitutes independent generalization.

## Structural posterior ranked fixed TRAIN predictions more effectively than raw gain

The completed whole-chromosome readout covered both strands of all six original TRAIN chromosomes for both adaptation branches, preserving their weights, MAP chains and the full 22,688-primary-reference denominator. Structural inclusion probabilities ranked primary matches more effectively than each model's raw chain gain in both species. Arabidopsis and rice AP increased from 0.7723 and 0.6816 to 0.9134 and 0.8832 for rank-four adaptation, and from 0.7542 and 0.6744 to 0.9147 and 0.8906 for full-rank adaptation. These are ranking measurements within each already-emitted prediction pool, not gene-detection AP over all references.

| Fixed TRAIN prediction pool | Raw gain AP | Structural posterior AP | Raw gain AUROC | Structural posterior AUROC |
| --- | ---: | ---: | ---: | ---: |
| Rank-four q/v | 0.7517 | 0.9078 | 0.5751 | 0.8130 |
| Full-rank q/v | 0.7350 | 0.9102 | 0.5521 | 0.8165 |

At the predeclared 90% retained-prediction budget, posterior ranking retained 16,046 primary matches among 20,092 rank-four predictions and 15,705 among 19,735 full-rank predictions. These were 689 and 634 more matches than gain ranking at the same respective counts, but 374 and 401 fewer than the unfiltered outputs. Recall against every primary reference therefore fell to 70.725% and 69.222%, respectively. No retained fraction was selected as an operating point, no filtered GFF was produced and no filtered CDS-base or background endpoint was inferred.

Structural inclusion was not calibrated biological correctness. Among predictions with structural probability at least 0.9, mean probabilities were 98.934% and 99.148%, whereas primary-reference match fractions were 85.879% and 85.034%; matches to any assessable annotated isoform were 92.314% and 91.798%. The corresponding pooled primary-reference Brier scores were 0.16345 and 0.17183, describing disagreement with those benchmark labels rather than validating the model's probability semantics. Ranking against any assessable isoform also improved, but reference-relative mismatches cannot all be classified as biologically false genes. The [completed readout and full accounting](../reports/M28-METHOD-RESTART/full_grammar_gene_posterior_train_plan.md) retain every species, fixed budget and probability bin. The separate missed-reference score assay was subsequently completed using four reconciled saved chromosomes and two unchanged-calculation recovery tasks, with the original timeout preserved in its [recovery record](../reports/M28-METHOD-RESTART/full_grammar_e_m_missed_reference_r2_recovery.md). This TRAIN-only exploration neither changes the frozen DEV families nor establishes independent generalization.

A separate completed readout of the unchanged equal-energy pool retained 16,527 primary matches among 22,066 predictions. Within this pool, posterior AP was 0.91017 versus 0.76157 for raw gain, and AUROC was 0.80535 versus 0.57264; both species showed the same ranking direction. The pooled posterior AUROC nevertheless remained below those of either parent on their respective prediction pools, and its primary-reference Brier score of 0.16651 lay between them. Different emitted geometries and positive prevalences prevent interpreting these differences as a paired probability-model superiority test. At the fixed 90% budget, posterior ranking retained 16,056 primary matches among 19,860 predictions, 603 more than gain ranking but 471 fewer than the unfiltered pool. No operating point was selected, no filtered GFF was generated and structural probabilities were not calibrated biological correctness. The [complete pool readout](../reports/M28-METHOD-RESTART/full_grammar_energy_pool_posterior_train.md) includes independent primary-label reconciliation of all chromosomes, species, fixed fractions and probability bins.

## Positive missed chain energies were often outscored by competing predictions

The completed unchanged-prediction score assay retained all 22,688 primary references and found identical single-chain grammar support for 22,479 in both adaptation branches. Of 6,268 rank-four misses, 3,896 (62.157%) had positive energy but were outscored by overlapping selected chains, 2,163 (34.509%) had negative energy against background and 209 (3.334%) were unsupported. Full-rank adaptation missed 6,582 references, with corresponding counts of 4,030 (61.228%), 2,343 (35.597%) and 209 (3.175%). Each legal replacement removed the entire gains of all overlapping same-strand selected chains; no improving replacement counterexample was found. Most missed references were therefore representable under the fixed single-chain grammar, but this score accounting does not establish why a different objective would recover them.

Both species retained this distinction between unsupported references and learned score competition. Negative start/stop contributions occurred even in 16,360 rank-four and 14,789 full-rank exact hits, so their sign alone could not identify the cause of a missed chain. The [full paired decomposition](../reports/M28-METHOD-RESTART/full_grammar_e_m_missed_reference_r2_recovery.md) preserves chromosome and species counts, retained single- and multisegment strata and all nine long-span references. Categories are descriptive reference-recall accounting, not subgroup precision or a causal role-weight experiment; predictions and the frozen DEV families remain unchanged.

## Separate completed native Carbon-A TRAIN comparison

The native Carbon-A pipeline was separately executed on the same six original TRAIN chromosomes and evaluated with the unchanged 22,688-primary-reference coordinate-and-strand ruler. Both predefined operating points were retained: confidence at least 0.1 or NA, and unfiltered predictions. No threshold was selected from these results, and neither method was added to the frozen DEV family.

| Fixed native operating point | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Worst species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: |
| Carbon-A, confidence at least 0.1 or NA | 78.958 | 75.259 | 78.172 | 77.108 | 75.259 |
| Carbon-A, unfiltered | 76.745 | 73.530 | 76.063 | 75.138 | 73.530 |

Filtered Carbon-A recovered 16,896 primary chains among 20,540 unique predictions, with pooled precision of 82.259% and recall of 74.471%. It exceeded ANNEVO, Tiberius and Helixer by 2.630, 5.516 and 14.580 percentage points pooled, and exceeded each in both species. Unfiltered Carbon-A recovered 17,023 chains among 22,072 predictions. Its pooled F1 exceeded ANNEVO by 0.522 points, but its Arabidopsis F1 was 0.086 points lower; a pooled advantage was not a uniform per-species advantage. The completed joint-encoding candidate remained 11.209 points below filtered Carbon-A pooled and 9.101 points below the unfiltered pipeline.

The fixed confidence filter removed 1,532 unique predictions and lost 127 primary hits, with no duplicate-chain removal. Pooled recall declined by 0.560 percentage points while F1 increased by 2.108 points. The other 1,405 removed chains were nonprimary-exact benchmark residuals, not necessarily biologically false genes. CDS-base F1 declined from 95.859% to 94.147%, while reference-relative background-span coverage declined from 4.545% to 1.282%. These endpoint changes describe the filter's tradeoff, not an independent causal effect or grounds to retune it.

The saved evaluator also distinguishes exact matches to other assessable annotated isoforms from primary hits. At the default filtered operating point, 1,365 of 3,644 nonprimary predictions exactly matched another assessable isoform; 2,279 matched neither the primary nor other assessable-isoform geometry. The corresponding unfiltered counts were 1,374 of 5,049 and 3,675, while the completed joint candidate had 1,267 of 6,966 and 5,699. Under the shorter native window policy, other-isoform matches numbered 1,361 filtered and 1,377 unfiltered. The default filter's 1,532 removed predictions comprised 127 primary hits, nine other-isoform matches and 1,396 matching neither reference category. These unique-prediction counts retain the unchanged primary endpoint; they are not per-gene recovery, explanations of individual primary losses, expression evidence or proof that unmatched structures are biologically false. The [reference-category decomposition](../reports/M28-METHOD-RESTART/carbon_a_original_train_result.md#非primary残差与已注释其他可评价-cds-链) preserves both species and operating points.

Reference-level replay found 14,336 chains recovered by both filtered Carbon-A and joint encoding, 2,560 Carbon-A-only hits, 590 joint-only hits and 5,202 missed by both. Multisegment references contributed 1,877 of the 1,970 net-hit difference (95.279%). Across one, two-to-five and at-least-six CDS segments, recall differences were 1.724, 9.133 and 13.158 percentage points, respectively. In rice, single-segment filtered Carbon-A had two fewer hits than joint encoding, whereas the other bins contributed net gains of 405 and 412. The unfiltered pipeline retained the same concentration: multisegment references accounted for 1,970 of its 2,097 net-hit gain.

Only nine primary references spanned more than 24,576 bases; joint encoding recovered none and either Carbon-A operating point recovered one. The other 22,679 references contributed 1,969 of the filtered pipeline's 1,970 net additional hits. Reference spans beyond the joint encoding window therefore cannot explain most of this panel's recovery deficit. This does not test the effect of longer context on shorter genes or identify a hard whole-contig grammar limit. These predeclared reference bins describe recall intersections, not subset prediction precision, F1, significance or an implemented ensemble.

A separate within-pipeline window-policy control reduced Carbon-A's neural window from 98,304 to 24,576 bases while retaining the released weights, all native postprocessing settings and the 6,144-base absolute overlap. Filtered pooled chain F1 declined from 78.172% to 77.562% (−0.609 percentage points), with 16,698 primary hits among 20,369 unique predictions; unfiltered F1 declined from 76.063% to 75.585% (−0.478 points), with 16,840 hits among 21,871 predictions. Both predefined operating points declined on all six chromosomes. At the shorter policy, filtered F1 was 78.453% in Arabidopsis and 74.292% in rice, with species-macro F1 of 76.372%; unfiltered values were 76.460%, 72.383% and 74.422%, respectively. The shorter native policy still exceeded the completed joint candidate by 10.600 and 8.622 points pooled. This intervention did not remove the main pipeline gap in the observed panel. Fixed absolute overlap changes relative overlap from 6.25% to 25% and stride from 92,160 to 18,432 bases, so the contrast combines window length, window placement and probability-averaging geometry rather than isolating context length. The default baseline is preserved, and this TRAIN-only [window-policy result](../reports/M28-METHOD-RESTART/carbon_a_context24_train_result.md) is not added to a frozen DEV family.

Carbon-A preserves the same annotation-pretrained main-weight source as the grammar model, but differs in context, native heads, learned transition and confidence networks, structural scoring and filtering. The contrast therefore does not isolate a decoder effect or establish independent generalization, clean zero-shot prediction or matched computational efficiency. The [completed native record](../reports/M28-METHOD-RESTART/carbon_a_original_train_result.md) preserves the results, unchanged predictions and transparent CPU recovery of one validator failure. The fresh-joint, uniform-coarse, codon-energy and low-rank encoder-adaptation TRAIN results are reported separately above; the completed full-rank continuation is reported separately above.

## Reference recovery, losses and additional-fit supervision

The [completed TRAIN-family figure](../reports/M28-METHOD-RESTART/figures/train_continuations/render_r1/train_continuation_complete_chain.svg) retains the initial model, all eight final continuations and three fixed native tools on shared per-species and pooled F1 axes. The [paired reference-hit figure](../reports/M28-METHOD-RESTART/figures/train_continuations/render_r1/train_continuation_reference_gain_loss.svg) separately displays gains and losses against the matched structured continuation, including negative branches. Their [captions and provenance](../reports/M28-METHOD-RESTART/figures/train_continuations/caption.md) preserve source jobs, counts, operating points and TRAIN-only interpretation; neither figure adds candidates to the frozen DEV comparison.

The completed reference replay retained all 22,688 primary TRAIN chains and all eight predeclared models, including the separate combined model. Exact-hit and unique-prediction counts reproduced the completed evaluations. Reference groups were fixed by the initial model and the three native tools, not redefined for each continuation. The 7,831-chain native gap comprised annotated primary references missed by the initial model but recovered by at least one native tool; 4,954 of these were recovered by all three tools. Recovering these chains is not discovery of new genes.

| Completed model | Recovered initial misses | Lost initial hits | Net primary hits | Recovered native gap / 7,831 | Recovered three-tool gap / 4,954 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Structured continuation | 3,366 | 636 | +2,730 | 3,235 | 2,483 |
| Local continuation | 2,119 | 2,542 | −423 | 1,995 | 1,468 |
| Mixed continuation | 2,938 | 1,692 | +1,246 | 2,763 | 2,002 |
| Native CDS prior | 4,906 | 508 | +4,398 | 4,753 | 3,571 |
| Fresh uniform draws | 5,159 | 985 | +4,174 | 4,970 | 3,703 |
| Joint 24,576-base encoding | 5,316 | 356 | +4,960 | 5,058 | 3,651 |
| Joint encoding plus native CDS prior | 4,629 | 906 | +3,723 | 4,495 | 3,364 |

Relative to the matched structured continuation rather than the lower-update initial model, joint encoding gained 2,956 references and lost 726, for a net increase of 2,230 hits. The prior gained 2,075 and lost 407; fresh sampling gained 2,421 and lost 977. Fresh sampling recovered more of the three-tool consensus gap than joint encoding, but had lower full-output chain F1 and more wholly background chains. A subset recovery count therefore does not identify the stronger complete-output predictor. These fixed outcome-defined groups provide descriptive recall accounting, not subset prediction precision, F1 or independent causal strata.

The predeclared four-cell comparison also did not support additive benefit from context and priors. On the pooled nonlinear F1 scale, joint encoding without priors exceeded the matched control by 11.217 percentage points, and priors without joint encoding by 8.594 points. Combining both exceeded the control by 6.706 points but fell below context alone by 4.512 points and priors alone by 1.888 points. The descriptive interaction, F11−F10−F01+F00, was −13.105 points pooled, −9.967 in Arabidopsis and −23.510 in rice. This single-seed TRAIN contrast does not establish statistical significance, mechanistic antagonism or equivalence of effective optimization.

For fresh versus original sampling, we stratified references by whether their whole span had directly matching fine labels in at least one sampled core during the additional fit. The original and fresh plans covered 2,097 and 4,505 references, respectively, with 595 in both. Fresh sampling's net gain over matched structured continuation was 1,444 hits: 39 in the both-covered stratum, 262 in fresh-only, 97 in original-only and 1,046 in neither. The neither-covered group contained 16,681 references and accounted for 72.438% of the net gain. This flag excludes common starting-fit exposure and does not imply zero encoder, partial-label, partition or upstream annotation exposure. The result does not support attributing all improvement to newly covered complete chains.

The complete-core proposal increased additional-fit coverage to 6,145 references, 1,640 more than fresh uniform sampling. Nevertheless, relative to fresh sampling it gained 1,531 and lost 1,774 exact-reference hits, for a net loss of 243: Arabidopsis gained 171, whereas rice lost 414. All four pooled coverage strata had negative net changes, including −90 in the 4,370-reference proposal-only stratum. Whole-output pooled F1 increased by 0.603 points, from 63.073% to 63.676%, while recall fell from 62.324% to 61.253% and precision rose from 63.840% to 66.299%; predictions fell from 22,149 to 20,961 and wholly background chains from 498 to 278. Increased complete-core exposure therefore did not produce uniformly higher recovery. All 209 grammar-unsupported references remained in the denominator. Importance correction matches the expected objective, not finite variance, clipped gradients or AdamW trajectories.

The [reference-recovery](../reports/M28-METHOD-RESTART/full_grammar_train_recovery_plan.md), [four-cell](../reports/M28-METHOD-RESTART/full_grammar_joint24_native_prior_plan.md) and [complete-core](../reports/M28-METHOD-RESTART/full_grammar_complete_core_draws_plan.md) records preserve the source jobs and fixed definitions. These supplemental analyses reuse completed TRAIN predictions, add no inference or training, and do not alter the frozen DEV candidates.

## Coarse observations retained primary region agreement under fine masks

The complete existing TRAIN-window pool was examined without changing labels, fits or predictions. Of 22,688 primary references, 22,452 were individually grammar-supported and fully contained in at least one existing core, giving 41,065 reference-core pairs. The 209 unsupported, 19 noncontained and eight overlong references remained unassessed rather than counted as agreements. Among the 6,220 references lacking any wholly matching fine core, 5,838 supplied 813,902 coarse-only row exposures: 378,329 coding-region and 435,573 intronic-region exposures. Every such exposure agreed with the primary reference region in both species. The remaining 382 references had no coarse-only observations in their saved containing cores. These repeated exposures are not unique genomic coverage or independent genes, and both labels and references derive from the same annotation.

For original coarse-only rows, clearing explicit fine-unknown intervals alone was sufficient to restore the primary fine role in 90,809 exposures (11.157%); removing alternatives alone was sufficient in 329,258 (40.454%). Another 389,158 (47.814%) required both changes, whereas 4,677 (0.575%) remained unresolved. No single-removal restoration disappeared under the joint change, and no newly known fine role under joint removal disagreed with the primary. These in-memory contrasts describe the target constructor, including changes to transcript-span occupancy; they neither identify incorrect isoforms nor demonstrate that removing masks would improve annotation. Fully unknown rows were not assigned to these coarse-only restoration classes. The [full assay record](../reports/M28-METHOD-RESTART/full_grammar_coarse_mask_census.md) retains every source category and the per-species counts. All actual training masks remained unchanged; this concordance result does not measure a partition gradient or complete-chain accuracy.

## Matched final update TRAIN scores explained residual omissions

A completed score replay compared structured continuation and joint encoding at the same final update 9,216, using their frozen six-chromosome TRAIN predictions. The assay was specified before the new DEV metrics were read, but not before the DEV jobs ran. It retained all 22,688 primary references, reproduced the unchanged prediction counts and selected-chain scores, and found no legal single-reference replacement that strictly improved either frozen selected set. References entered only post hoc path scoring, not the encoder, head or inference. The [completed assay record](../reports/M28-METHOD-RESTART/full_grammar_matched_energy_diagnostic.md#完整终态与评分结果) contains the native replay and arithmetic checks.

Joint encoding missed 7,762 primary chains: 2,486 had negative reference energy relative to background, 5,067 had positive energy but were outscored by the entire set of selected same-strand conflicts, and 209 were outside the single-chain grammar. Thus 97.307% of its omissions were individually grammar-supported; this does not mean that all supported references could coexist in one compatible set. Relative to structured continuation, the negative-energy category decreased by 1,870, but the competition category decreased by only 360 pooled and increased by 125 in rice.

| TRAIN scope | Structured negative energy | Joint negative energy | Structured outscored | Joint outscored | Individually unsupported, both |
| --- | ---: | ---: | ---: | ---: | ---: |
| Arabidopsis | 2,629 | 1,667 | 4,165 | 3,680 | 108 |
| Rice | 1,727 | 819 | 1,262 | 1,387 | 101 |
| Pooled | 4,356 | 2,486 | 5,427 | 5,067 | 209 |

The paired migration distinguished recovery from reclassification. Of the 2,956 newly matched chains, 2,123 (71.820%) had previously been positive-energy references outscored by selected conflicts, and 833 had negative energy. A further 1,085 references moved from negative energy to positive-but-outscored without becoming exact hits. Joint encoding also lost 726 hits, yielding the unchanged net gain of 2,230. A decrease in the negative-energy category alone therefore does not measure recovered complete chains.

In the 2,956 gained-chain stratum, the structured model's mean replacement margin was −18.4933 score units: coding body, intron body, splice-coding boundary and start/stop contributions were +5.1171, −29.0444, +2.9378 and +2.4963. In the 6,827 supported chains missed by both models, the joint model's mean margin was −47.4246, with corresponding contributions −7.2952, −34.7842, −2.4944 and −2.8508. Its intron contribution was −20.1887 in Arabidopsis and −67.7967 in rice; the rice splice-boundary contribution was nevertheless positive (+3.4700). Intron terms were also the largest negative mean contribution for the 726 chains lost by joint encoding (−11.0849), despite a positive splice-boundary contribution (+5.5680). These are within-model score budgets, not comparable probability scales across heads.

The accounting identifies an intron-associated score deficit for these fixed outcome groups, not a causal failure of intron modeling or a reason to reweight inference scores. Contributions are unnormalized sums over differently sized reference and competing chains, not per-base error rates or unique genomic coverage. A selected exact chain has zero replacement margin by identity, which cannot establish the mechanism of its recovery. Outcome-defined groups, a single fit and reused TRAIN data preclude causal or independent-generalization claims. Neither the grammar support counts nor this score arithmetic erase the remaining performance gap to ANNEVO and Tiberius.

## Coarse observation signal varied across fixed TRAIN views

We held the joint-encoding update-9,216 weights and native logits fixed and compared the informative coarse channel with a uniform channel and the same fine observations alone. All three calculations used each draw's original fine-plus-coarse denominator. The first coarse-bearing draw per species was fixed before this diagnostic; after those results were known, a separate follow-up selected each species' maximum coarse count from the unchanged fresh schedule, with earliest draw breaking ties. The latter selection uses label counts, not native loss, gradients, predictions or DEV. Both read-only jobs completed without parameter updates or prediction changes.

| Species | Selection | Coarse-only rows | Original denominator | Normalized informative-minus-uniform logit-gradient L2 |
| --- | --- | ---: | ---: | ---: |
| Arabidopsis | First coarse-bearing draw | 1 | 24,480 | 0 |
| Rice | First coarse-bearing draw | 3 | 24,357 | 0 |
| Arabidopsis | Maximum coarse count | 7,127 | 24,576 | 0.001013543397 |
| Rice | Maximum coarse count | 17,474 | 24,576 | 0.005446409732 |

Uniform minus fine-only loss matched the expected coarse-row entropy constant, with absolute residuals at most 1.46×10⁻¹¹ and identical returned logit derivatives in all four views. The two initial informative losses also differed from fine-only by the all-correct coarse constant, whereas the maximum-count views had measurable additional derivatives. Thus coarse observations can supply a conditional learning signal at this fixed head; the two zero-gradient examples do not establish universal redundancy. These selected points are not a population sample, parameter-gradient analysis, training ablation or accuracy measurement. Partitions used FP64 and derivatives returned at native FP32 logit precision. The separately fixed matched continuation has completed: informative minus uniform chain F1 was -0.908 points in Arabidopsis, +1.537 in rice and -0.342 pooled. Local nonzero derivatives therefore did not imply a uniform performance benefit. The [diagnostic record](../reports/M28-METHOD-RESTART/full_grammar_coarse_information_probe.md) retains counts, loss identities, role-group derivatives and selection timing, and the [continuation plan](../reports/M28-METHOD-RESTART/full_grammar_uniform_coarse_continuation_plan.md) retains the common already-trained starting point and its limitation.

## Three continuation objectives on the fixed DEV panel

We evaluated all three predeclared final continuation weights on the two fixed DEV chromosomes, totaling 53,396,251 bases. Candidate, native-comparator and original FASTA records were directly equal before reference access. The unchanged common ruler retained 5,437 Arabidopsis and 2,291 rice primary chains, including unsupported reference strata. The original final-4,608 anchor and six native outputs were reused without new native runs. Each inference task and its corresponding evaluation completed successfully; pooled counts, base metrics, background totals, species macro and worst-species scores reconciled, and the reused native measurements were identical across the four reports.

| Model or fixed comparator | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Wholly background chains |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial structured model, 4,608 updates | 50.844 | 30.096 | 44.608 | 40.470 | 64 |
| Structured continuation, 9,216 updates | 62.923 | 36.931 | 54.751 | 49.927 | 79 |
| Local continuation, 9,216 updates | 44.485 | 36.867 | 42.181 | 40.676 | 173 |
| Mixed continuation, 9,216 updates | 51.801 | 42.452 | 48.943 | 47.127 | 245 |
| ANNEVO | 76.265 | 74.432 | 75.703 | 75.348 | 146 |
| Tiberius | 74.177 | 66.582 | 71.805 | 70.380 | 165 |
| Helixer | 67.456 | 53.524 | 62.851 | 60.490 | 563 |

Structured continuation recovered 4,267 primary chains among 7,859 predictions: pooled precision was 54.294%, recall was 55.215% and F1 was 54.751%, 10.142 percentage points above the initial anchor. It exceeded local continuation by 12.570 points pooled despite 1.051-point lower CDS-base F1. Mixed continuation exceeded structured continuation in rice by 5.521 points, but was 11.122 points lower in Arabidopsis and 5.808 points lower pooled. Its pooled CDS-base F1 was 93.895%, compared with 92.756% for structured continuation, and it produced 166 more wholly background chains. Rice remained the worst species for all seven methods.

All three continuation objectives remained below all three native comparators in both species. The highest pooled continuation score was still 8.100 percentage points below Helixer, 17.054 below Tiberius and 20.952 below ANNEVO. The seven-method table is the completed objective subcomparison; the complete eight-method family and separate ten-method extension are reported below. This historically reused same-species development panel cannot establish independent generalization, clean zero-shot performance or publication competitiveness, and its readings are not used to select checkpoints or change subsequent training. The [frozen objective DEV record](../reports/M28-METHOD-RESTART/full_grammar_objective_contrast_dev_plan.md#三目标-dev-终态与结果) retains the counts, terminal receipts and protocol.

## Complete frozen DEV family and predeclared extension

The native-prior, fresh-sampling and joint-encoding branches completed both whole-chromosome DEV inference tasks and their common-ruler evaluations. The original eight-method aggregation retained the initial model, three objective continuations, native-prior continuation and three fixed native tools. A separate predeclared extension appended fresh sampling and joint encoding without changing any of those eight measurements. The combined and complete-core branches were not added. All seven M28 source comparisons retained the same 7,728 primary references, unsupported and long-span strata, and directly identical native measurements. Count-based chain and CDS-base precision, recall and F1, background totals, species macro and worst-species F1, phase and structural partitions reconciled independently.

The reviewed [complete-chain F1 figure](../reports/M28-METHOD-RESTART/figures/frozen_dev_family/render_r1/dev_family_complete_chain.svg) and [fixed-operating-point precision/recall figure](../reports/M28-METHOD-RESTART/figures/frozen_dev_family/render_r1/dev_family_precision_recall.svg) show the complete original eight-method family and both predeclared extensions, retaining the fixed order and common 0–100% axes. Their [captions and export record](../reports/M28-METHOD-RESTART/figures/frozen_dev_family/caption.md) preserve source counts, completion evidence and the development-panel limitations. The precision/recall bars distinguish an F1 advantage from a recall advantage; they are not a threshold sweep.

| Fixed model or comparator | Arabidopsis chain F1 % | Rice chain F1 % | Pooled chain F1 % | Species macro F1 % | Wholly background chains |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial structured model, 4,608 updates | 50.844 | 30.096 | 44.608 | 40.470 | 64 |
| Structured continuation, 9,216 updates | 62.923 | 36.931 | 54.751 | 49.927 | 79 |
| Local continuation, 9,216 updates | 44.485 | 36.867 | 42.181 | 40.676 | 173 |
| Mixed continuation, 9,216 updates | 51.801 | 42.452 | 48.943 | 47.127 | 245 |
| Native CDS prior, 9,216 updates | 69.297 | 51.516 | 63.776 | 60.407 | 121 |
| ANNEVO | 76.265 | 74.432 | 75.703 | 75.348 | 146 |
| Tiberius | 74.177 | 66.582 | 71.805 | 70.380 | 165 |
| Helixer | 67.456 | 53.524 | 62.851 | 60.490 | 563 |
| Fresh uniform draws, predeclared extension | 67.438 | 51.909 | 62.685 | 59.673 | 188 |
| Joint 24,576-base encoding, predeclared extension | 70.499 | 57.502 | 66.623 | 64.000 | 121 |

Rice was the worst species for every method, so its column also gives worst-species F1. Native priors and fresh sampling exceeded matched structured continuation by 9.025 and 7.934 percentage points pooled, respectively. The prior exceeded Helixer by 0.925 points pooled but was 2.008 points lower in rice and 0.084 points lower on the species macro endpoint. Fresh sampling remained below Helixer in both species. An aggregate advantage therefore did not establish a uniform advantage for either branch.

Joint encoding recovered 5,070 primary chains among 7,492 unique predictions, with pooled precision of 67.672%, recall of 65.606% and F1 of 66.623%. Relative to matched structured continuation, it added 803 net primary hits and improved F1 by 7.576 points in Arabidopsis, 20.570 in rice and 11.872 pooled. It exceeded Helixer's chain F1 by 3.043 and 3.977 points in the two species and 3.772 pooled. However, Helixer recovered 5,112 references among 8,539 predictions: joint encoding had 42 fewer hits and 1,047 fewer predictions, with 7.806-point higher pooled precision but 0.543-point lower recall. In Arabidopsis it had 92 more hits; in rice it had 134 fewer. The F1 advantage is not a claim of higher recall in both species.

Joint encoding remained below Tiberius by 3.678 points in Arabidopsis, 9.081 in rice and 5.182 pooled, and below ANNEVO by 5.766, 16.930 and 9.080 points, respectively. Its pooled CDS-base F1 was 95.096%, compared with 93.672% for native priors, 93.994% for fresh sampling and 92.031% for Helixer; base coverage does not replace the complete-chain endpoint. Its 121 wholly background chains were 442 fewer than Helixer but 42 more than structured continuation. All methods retained the three references spanning more than a native read window, and all seven M28 weights missed those three; that small descriptive stratum does not establish long-gene capability.

These were the predeclared final weights, not candidates selected by DEV performance. The two historically reused same-species DEV chromosomes, single fit seed and method-dependent upstream annotation exposure remain important limitations. Neither the superiority to one comparator nor the transfer of a TRAIN improvement to this development panel establishes independent generalization, a controlled causal mechanism, novelty or publication readiness. The [original eight-method record](../reports/M28-METHOD-RESTART/full_grammar_native_cds_prior_dev_plan.md#完整八方法-dev-终态与结果) and [separate ten-method extension](../reports/M28-METHOD-RESTART/full_grammar_sampling_context_dev_plan.md#十方法扩展-dev-终态与结果) retain all source jobs, counts and unchanged protocol.


## Three fixed fresh-joint candidates on the reused DEV panel

All three fixed final-update-9,216 candidates completed raw-FASTA inference and common-ruler evaluation on the same two historical DEV contigs. The separate six-method family retained fresh-joint, its uniform-coarse control and its codon-energy extension regardless of TRAIN ranking, alongside the unchanged native outputs. Direct caller-FASTA equality was established before reference access, and all 7,728 primary references, including unsupported and long-span strata, remained in evaluation. Earlier eight- and ten-method results were not changed.

| Method | Arabidopsis F1 % | Rice F1 % | Pooled primary hits / unique predictions | Pooled F1 % | Species macro F1 % | Worst-species F1 % |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Fresh-joint | 72.858 | 63.725 | 5,353 / 7,543 | 70.107 | 68.292 | 63.725 |
| Uniform coarse | 73.318 | 61.641 | 5,355 / 7,617 | 69.795 | 67.480 | 61.641 |
| Codon energy | 69.038 | 54.330 | 5,041 / 7,884 | 64.579 | 61.684 | 54.330 |
| ANNEVO | 76.265 | 74.432 | 5,644 / 7,183 | 75.703 | 75.348 | 74.432 |
| Tiberius | 74.177 | 66.582 | 5,441 / 7,427 | 71.805 | 70.380 | 66.582 |
| Helixer | 67.456 | 53.524 | 5,112 / 8,539 | 62.851 | 60.490 | 53.524 |

Fresh-joint recovered 5,353 primary chains among 7,543 predictions, with pooled precision of 70.966%, recall of 69.268% and F1 of 70.107%. Relative to Helixer, it added 241 net primary hits and made 996 fewer predictions; precision, recall and F1 were higher in both species. Its pooled F1 advantage was 7.256 points. It nevertheless remained 1.698 points below Tiberius and 5.596 below ANNEVO, with deficits in both species. All three grammar candidates had higher F1 than Helixer and lower F1 than the two stronger native comparators; this does not establish the competitive goal.

The informative coarse channel exceeded its uniform control by 0.312 points pooled, but was 0.459 points lower in Arabidopsis and 2.084 higher in rice. It recovered two fewer primary chains pooled while making 74 fewer predictions. The species tradeoff persisted, and the pooled contrast had the opposite sign from the completed TRAIN contrast; neither view demonstrates a universal or full-training necessity of coarse information. Codon energies lowered F1 relative to fresh-joint by 3.820 points in Arabidopsis, 9.395 in rice and 5.528 pooled, with 312 fewer primary hits and 341 more predictions pooled. This single-fit intervention remained a negative result in both the original TRAIN range and reused DEV panel.

Pooled CDS-base F1 was 95.065%, 94.960% and 94.859% for fresh-joint, uniform and codon, respectively, while wholly-background chains numbered 145, 154 and 196. Fresh-joint had higher base F1 than Tiberius and ANNEVO despite lower complete-chain F1, so base coverage cannot substitute for the primary endpoint. It had fewer wholly-background chains than Tiberius but higher background-span coverage, 0.921% versus 0.721%. All three candidates missed the three retained references spanning more than a native read window; this small stratum establishes no long-gene capability.

The [completed fixed-family record](../reports/M28-METHOD-RESTART/full_grammar_joint_controls_dev_result.md) retains every candidate, full counts, base/background endpoints, other-isoform and phase accounting, and negative comparisons. These are single-fit measurements on historically reused same-species DEV contigs, not independent generalization or clean zero-shot evidence. No DEV reading selects checkpoints, parameters, subsequent training or operating points; both completed adaptation branches remain outside this family.

## Joint structural compatibility did not explain the large accuracy gap

We assessed whether the grammar could represent the reference chains jointly, rather than only one chain at a time. Of the 22,688 primary references, 22,479 satisfied the implemented DNA grammar individually. For each chromosome and strand, we then found the maximum number of supported reference CDS spans that could coexist without overlap. Touching half-open spans were compatible, and opposite strands were evaluated independently. The maximum compatible set contained 22,460 references: 209 references were individually unsupported, and joint compatibility imposed a further loss of 19 supported references.

| TRAIN scope | All references N | Maximum supported compatible references M | Minimum forced omissions | Reference oracle F1 % |
| --- | ---: | ---: | ---: | ---: |
| Arabidopsis | 17,739 | 17,613 | 126 | 99.6436 |
| Rice | 4,949 | 4,847 | 102 | 98.9588 |
| Pooled | 22,688 | 22,460 | 228 | 99.4950 |

If an oracle outputs only these M reference chains and no additional predictions, its precision is one, recall is M/N and F1 is 2M/(N+M). The pooled value of 99.4950% is therefore a reference-informed structural ceiling within the declared grammar. It is not raw-FASTA model accuracy, a guarantee that the neural scores can realize the oracle, or evidence of independent generalization. A reference omitted from one maximum set may belong to another equally large set.

All 16,470 primary references correctly recovered by ANNEVO and all 16,095 recovered by Tiberius were individually grammar-supported and jointly compatible within each tool's own hit set. This statement does not combine the two tools' hit sets. Helixer recovered 15,091 primary references, of which 15,072 were grammar-supported and jointly compatible; its remaining 19 hits were outside the grammar and stayed in the original evaluation. These results show that the grammar and same-strand nonoverlap constraint alone cannot explain the initial model's large deficit relative to ANNEVO and Tiberius. They do not identify whether scoring, supervision, optimization or their interaction is responsible. The [joint compatibility diagnostic](../reports/M28-METHOD-RESTART/full_grammar_train_energy_diagnostic.md#原train参考集合的联合结构上限) records the exact definitions, completed oracle calculation and retained denominators.

## Scope of the evidence

The [Discussion draft](M28_full_grammar_discussion_working.md) relates these complete-chain measurements and score diagnostics to the contribution under test, without treating them as a causal mechanism or independent validation.

These are single-fit diagnostics on the initial training range and historically reused same-species DEV chromosomes, not independent tests or a matched training-budget comparison with the native tools. Upstream pretraining and annotation exposure are model-dependent; freezing an annotation-pretrained backbone does not establish exposure-free prediction. The [exposure record](../reports/M28-METHOD-RESTART/pretrained_exposure.md) preserves these distinctions. The scope is complete protein-coding CDS chains, not UTR annotation or recovery of every alternative isoform.

The completed inference runs also do not establish a computational-efficiency advantage. CPU allocations, software environments, concurrent node workloads and timing boundaries differed between the initial model and the native wrappers; the [execution record](../reports/M28-METHOD-RESTART/native_train_baseline_result.md#完成运行的资源占用) reports the observed allocations without treating them as a fair speed benchmark. Competitive complete-chain performance, a controlled explanation of any improvement and independently supported broader applicability remain to be established by the corresponding experiments.
