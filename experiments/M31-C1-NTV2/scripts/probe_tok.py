import time, torch
from transformers import AutoTokenizer, AutoModelForMaskedLM
M = "InstaDeepAI/nucleotide-transformer-v2-500m-multi-species"
tok = AutoTokenizer.from_pretrained(M, trust_remote_code=True)
print('vocab', len(tok.get_vocab()), [k for k, v in tok.get_vocab().items() if len(k) != 6])
for s in ['ACGTAC' * 3, 'ACGTAC' * 2 + 'ACGNAC', 'N' * 12, 'ACGTACG']:
    print(s, tok(s, add_special_tokens=True)['input_ids'], tok(s, add_special_tokens=False)['input_ids'])
