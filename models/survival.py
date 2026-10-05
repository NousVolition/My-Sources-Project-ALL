"""Series survival checks from Survives, Survives2, Ways to Failure, Translated Failures."""


def simulate_balance_metaphor(body, mind, symmetry, tension, utility):
    return body * mind * symmetry * tension * utility


def check_survival(B, M, Y, T, U):
    if B and M and Y and T and U:
        return "System Active: Metaphor survives! (1)"
    return "System Inactive: Metaphor fails. (0)"


def check_system_state(E, P, S, T, U):
    B = E * P * S * T * U
    if B == 1:
        return B, "Nominal: Metaphor survives."
    failures = []
    if E == 0:
        failures.append("Disembodied (E=0)")
    if P == 0:
        failures.append("Literalized (P=0)")
    if S == 0:
        failures.append("Asymmetric (S=0)")
    if T == 0:
        failures.append("Static (T=0)")
    if U == 0:
        failures.append("Obsolete (U=0)")
    return B, f"Failure due to: {', '.join(failures)}"


def check_aerospace_state(B, M, Y, T, U):
    S = B * M * Y * T * U
    if S == 1:
        return S, "Nominal Operation: Stable Flight Profile Verified."
    failures = []
    if B == 0:
        failures.append("Airframe Compromise (B=0)")
    if M == 0:
        failures.append("Avionics Dark (M=0)")
    if Y == 0:
        failures.append("Asymmetry Condition (Y=0)")
    if T == 0:
        failures.append("Thrust Stall (T=0)")
    if U == 0:
        failures.append("Payload Inefficient (U=0)")
    return S, f"System Failure due to: {', '.join(failures)}"
