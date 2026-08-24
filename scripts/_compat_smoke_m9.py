"""M9-CK1 compatibility smoke: can transformers 5.11.0 load NT-v2-500m via trust_remote_code,
expose a trainable .esm backbone (1024-dim hidden), and forward a tiny batch? Also probe the
6-mer tokenizer (token->base ratio) for the per-base alignment. CPU forward of a short seq.
Run: python scripts/_compat_smoke_m9.py    (prints COMPAT_OK / COMPAT_FAIL:<reason>)
"""
import sys, traceback
MODEL = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
try:
    import torch, transformers
    from transformers import AutoModelForMaskedLM, AutoTokenizer
    print("transformers", transformers.__version__, "torch", torch.__version__)
    tok = AutoTokenizer.from_pretrained(MODEL, trust_remote_code=True)
    # tokenizer probe: encode a 60bp ACGT seq, see how many tokens (6-mer -> ~10 tokens + specials)
    seq = "ACGT" * 15  # 60 bp
    enc = tok(seq, return_tensors="pt")
    ntok = enc["input_ids"].shape[1]
    print(f"tokenizer: 60bp -> {ntok} tokens (6-mer expect ~10-12 incl specials); ids[:8]={enc['input_ids'][0,:8].tolist()}")
    model = AutoModelForMaskedLM.from_pretrained(MODEL, trust_remote_code=True, output_hidden_states=True)
    print("loaded AutoModelForMaskedLM; type", type(model).__name__)
    esm = getattr(model, "esm", None)
    print("has .esm backbone:", esm is not None, "| type", type(esm).__name__ if esm else None)
    if esm is not None:
        nlayers = len(esm.encoder.layer)
        print(f".esm.encoder.layer: {nlayers} layers")
    model.eval()
    with torch.no_grad():
        out = model(**enc)
    hs = out.hidden_states[-1] if getattr(out, "hidden_states", None) is not None else None
    print("forward OK; last hidden_state shape:", tuple(hs.shape) if hs is not None else "NO hidden_states")
    # gradient_checkpointing available?
    gc = hasattr(model, "gradient_checkpointing_enable")
    print("gradient_checkpointing_enable available:", gc)
    if hs is not None and hs.shape[-1] == 1024 and esm is not None:
        print("COMPAT_OK")
    else:
        print("COMPAT_FAIL: hidden dim or .esm unexpected")
except Exception as e:
    print("COMPAT_FAIL:", repr(e))
    traceback.print_exc()
