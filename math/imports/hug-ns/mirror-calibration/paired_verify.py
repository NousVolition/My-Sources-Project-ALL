"""Complete paired symmetry and signed calibration tests; reuse completed runs."""
from mirror_calibration import main
from analyze_mirror_calibration import rhs_checks, analyze, plots

if __name__ == '__main__':
    main()
    rhs_checks()
    result = analyze()
    plots(result)
    print('Results: outputs/mirror-calibration.html')
