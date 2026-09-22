"""Format count-derived estimates without rounding a rounded table estimate."""
from decimal import Decimal, ROUND_HALF_UP, localcontext


def format_mcc_from_counts(row, places=3):
    """Compute MCC at high precision from TP/FP/TN/FN, then round once."""
    tp, fp, tn, fn = (Decimal(str(row[key])) for key in ("TP", "FP", "TN", "FN"))
    with localcontext() as context:
        context.prec = 40
        product = (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)
        value = (tp * tn - fp * fn) / product.sqrt() if product else Decimal(0)
        return str(value.quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP))
