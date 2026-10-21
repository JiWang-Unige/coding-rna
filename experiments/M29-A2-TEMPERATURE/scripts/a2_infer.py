"""A2 inference = A1 inference with the bias-free A2 head class (same offset, same decoder). Usage identical to a1_infer.py."""
import a2_common as C2     # noqa: F401
import a1_model
import a2_model
a1_model.new_model = a2_model.new_model
import a1_infer

if __name__ == '__main__':
    a1_infer.main()
