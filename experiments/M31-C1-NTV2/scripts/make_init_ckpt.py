"""Random-init head checkpoint (engineering test of the inference/scoring path only)."""
import sys, torch
import c1_model as MD
layers = (21, 27)
m = MD.new_head(layers)
torch.save(dict(model=m.state_dict(), step=0, config=dict(layers=list(layers), layer_norm_per_layer=True)), sys.argv[1])
