"""Random-init head checkpoint (engineering test of the inference/scoring path only)."""
import sys, torch
import c1_model as MD
m = MD.new_head()
torch.save(dict(model=m.state_dict(), step=0, config=dict(layer=-1)), sys.argv[1])
