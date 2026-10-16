"""Read out the probability mass assigned to a specified legal event support."""
import torch


def conditional_support(log_weights, target_index):
    if not len(log_weights) or not 0 <= target_index < len(log_weights):
        raise ValueError('A nonempty support containing the actual target is required')
    values = log_weights.double()
    log_mass = torch.logsumexp(values, dim=0)
    return {'log_mass': log_mass, 'conditional_target_logp': values[target_index]-log_mass}
