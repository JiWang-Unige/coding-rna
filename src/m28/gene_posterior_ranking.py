"""Reference-relative diagnostics of already selected CDS chains; no fitting."""
import math


def oriented_chain(parts, length, strand):
    parts = tuple((int(a), int(b)) for a, b in parts)
    if strand == '+':
        return parts
    if strand == '-':
        return tuple((length-b, length-a) for a, b in reversed(parts))
    raise ValueError('One genomic strand required')


def ranking_summary(rows, score, label, reference_count=None):
    """Group exact score ties. AP/AUC classify the selected prediction pool only.

    Primary-reference recall has its separate full retained denominator; AP is
    NOT annotation detection AP and excludes all unpredicted reference genes.
    """
    if not rows or any(not math.isfinite(r[score]) or type(r[label]) is not bool for r in rows):
        raise ValueError('Nonempty finite scores and explicit binary reference-relative labels required')
    ordered = sorted(rows, key=lambda r: -r[score])
    n, positive = len(rows), sum(r[label] for r in rows)
    negative = n-positive
    seen, tp, ap, concordant = 0, 0, 0., 0.
    curve, cursor = [], 0
    while cursor < n:
        end = cursor+1
        while end < n and ordered[end][score] == ordered[cursor][score]:
            end += 1
        group = ordered[cursor:end]
        gp = sum(r[label] for r in group); gn = len(group)-gp
        seen += len(group); tp += gp
        ap += gp*tp/seen
        # Positive beats every negative below this group, ties count one half.
        concordant += gp*(negative-(seen-tp)+.5*gn)
        curve.append(dict(score=ordered[cursor][score], retained_predictions=seen,
            true_reference_matches=tp, precision=tp/seen,
            recall_within_selected_positive_pool=tp/positive if positive else None,
            recall_of_all_primary_references=tp/reference_count if reference_count else None))
        cursor = end
    retained = {}
    for fraction in (.25, .5, .75, .9, 1.):
        row = next(r for r in curve if r['retained_predictions'] >= math.ceil(n*fraction))
        retained[str(fraction)] = dict(row, requested_fraction=fraction,
            actual_fraction=row['retained_predictions']/n, ties_kept_together=True)
    return dict(predictions=n, positive=positive, negative=negative,
        prediction_pool_average_precision=ap/positive if positive else None,
        prediction_pool_AUROC=concordant/(positive*negative) if positive and negative else None,
        retained_fraction_diagnostics=retained, grouped_precision_recall_curve=curve,
        unpredicted_reference_genes_in_AP_or_AUROC=False, threshold_selected=False)


def probability_diagnostics(rows, label):
    """Brier and ten fixed equal-width bins against a reference-relative event.

    The probability models structural inclusion, not this correctness event.
    This diagnostic can expose mismatch; it does not establish biological
    calibration, and nonmatching predictions are not proven false genes.
    """
    bins = [dict(lower=i/10, upper=(i+1)/10, count=0, probability_sum=0., matches=0) for i in range(10)]
    brier = 0.
    for row in rows:
        p = math.exp(row['log_probability'])
        if not 0 <= p <= 1 or type(row[label]) is not bool:
            raise ValueError('Valid structural probability and binary diagnostic event required')
        y = row[label]
        brier += (p-y)**2
        bucket = bins[min(9, int(p*10))]
        bucket['count'] += 1; bucket['probability_sum'] += p; bucket['matches'] += y
    for bucket in bins:
        n = bucket['count']
        bucket['mean_structural_probability'] = bucket.pop('probability_sum')/n if n else None
        bucket['reference_exact_fraction'] = bucket['matches']/n if n else None
    return dict(reference_relative_Brier=brier/len(rows), fixed_probability_bins=bins,
        probability_target='structural_inclusion_not_reference_exact_correctness',
        empirical_biological_calibration_claimed=False, fitted_temperature=False)
