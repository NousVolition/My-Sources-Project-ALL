"""Prespecified settings for the three newly requested configurations."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GAMMA=12.138737596974998
GRIDS=(48,64,80,112,160)
CASES=('aligned','compressive','exodus')
METHODS=('fourier','fd4')
NU=.001
END=.40
OUTPUT_DT=.02
CFL_LIMIT=.75

def source_hashes():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [ROOT/'initial_design.py',ROOT/'solver.py',ROOT/'protocol.py',ROOT/'run_study.py']}

def jobs(case=None):
    for c in CASES if case is None else (case,):
        for n in GRIDS:
            for half in ((False,True) if n>=80 else (False,)):
                for method in METHODS:
                    steps_per_output=round(25*n/80)*(2 if half else 1)
                    yield {'id':f'{c}-{method}-n{n}-'+('half' if half else 'base'),
                           'case':c,'method':method,'n':n,'half':half,'nu':NU,'end':END,
                           'output_dt':OUTPUT_DT,'steps_per_output':steps_per_output,
                           'dt':OUTPUT_DT/steps_per_output,'gamma':GAMMA}

def protocol():
    return {
      'title':'Aligned, compressive, and departing vortex surroundings',
      'scope':'New initial construction. These runs do not continue the original ring, hug demo, or imported matched-stretch fields.',
      'equation':'Periodic incompressible Navier-Stokes; external force identically zero; density one.',
      'domain_side':6.,'viscosity':NU,'end_time':END,'output_interval':OUTPUT_DT,
      'design_version':2,
      'chosen_parameters':{'central_vorticity_peak':80.,'central_gaussian_scale_radius':.70,
         'central_axial_scale':1.,'outer_scale_radius':.70,'outer_gamma':GAMMA,'outer_y_centres':[-1.5,1.5],
         'aligned_and_compressive_outer_z_centres':[-1.5,1.5],'exodus_outer_z_centres':[-.9,.9]},
      'construction':[
        'Smooth analytic periodic cosine-exponential profiles; central localized tube along z, four surrounding tubes along x.',
        'The central velocity is curl(0,0,psi), where psi = 40 * 0.70^2 * bump(x,y;0.70) * axial_bump(z;1.0). It has smooth return and transverse vorticity; no tube ends or discontinuities.',
        'Both formulas use scale 0.70, but the central vorticity is a second derivative of its bump and is narrower than the outer vorticity bumps. Resolution is assessed from measured half-peak chords, not from the scale label.',
        'The outer velocity uses periodic Biot-Savart. All components have zero mean and are divergence-free. No imposed strain or other force remains after initialization.',
        'Gamma was calibrated once at 80 cubed for initial axial background strain +4 at the origin; it is then fixed for all grids and cases.',
        'Aligned and compressive differ only by sign of the background. Exodus also changes the outer z spacing to create outward initial motion.',
        'The compressive variant imposes initial compression, not enforced tilting, bending, dispersion, or later decay.',
        'Gaussian tails overlap; well-separated means separated half-peak cores initially, not disjoint support.',
        'These selected parameters define an adversarial experiment; maximality among all admissible flows is not established.'
      ],
      'rejected_pilot':'Version 1 used an axially uniform central tube; reversing its quadrupolar background was exactly a half-box translation. That pilot was stopped and preserved under design-v1-rejected. Version 2 localizes the central tube axially, and its initial common-location stretching response is checked with the full RHS.',
      'discretization':{
        'fourier':'Rotational transport, spectral derivatives and Laplacian, Fourier pressure projection.',
        'fd4':'Physical-space fourth-order finite-difference skew transport and Laplacian; native discrete divergence and gradient. FFT diagonalizes its discrete pressure Poisson operator.',
        'shared':'Strict rectangular |mode| < N/3 cutoff and explicit SSP RK3 time integration. No peak cap, body force or feedback.',
        'independence_limit':'FD4 independently discretizes momentum transport and viscosity, but shares FFT pressure infrastructure, filtering, RK3 and diagnostics. It is not an entirely FFT-free or independently authored code.',
        'fd_initial_projection':'The common analytic velocity is projected to the FD divergence-free subspace. The relative L2 change is recorded on each grid.'},
      'measurements':{
        'W':'Maximum magnitude of Fourier curl of interpolated velocity over sampled grid points, shared diagnostic across methods; also record native FD curl.',
        'I':'Integrate maximum spin at every RK stage with RK3 weights. This numerical integral is not an upper bound on the continuum supremum.',
        'central_roi':'Fixed radius-1 cylinder about original z axis; complements global peak and passive central markers. It is not automatic vortex identity.',
        'width':'12 transverse half-peak chords through a sampled peak, in physical domain units and grid cells. Subcell interpolation does not add resolution.',
        'markers':'Passive Heun-tracked material points; periodic nearest distances between initially central and outer centrelines. Viscosity can separate vorticity from material markers.',
        'budgets':'Volume-integrated energy and native enstrophy, stage-integrated production/loss and balance residuals, direct physical stretching integral at outputs.',
        'thresholds':[50,100,200,400],
        'units':'Model length and time units; no physical SI calibration.'},
      'numerical_gates':{'max_advective_stage_cfl':CFL_LIMIT,'scaled_native_divergence':1e-10,
        'absolute_relative_energy_balance':.005,'minimum_initial_core_chord_cells':8.,
        'stop':'Nonfinite state, CFL breach, failed incompressibility, or energy balance above tolerance stops the matrix and preserves the field. These are numerical validity gates, not physical opening rules.',
        'warnings':'Thin measured chords below 6 cells or high-band enstrophy above 1 percent flag unresolved structure but do not replace refinement comparisons.',
        'comparison':'Report all differences before describing agreement; 5 percent peak/integral and width differences are screening tolerances, not proof of asymptotic convergence.'},
      'case_order':list(CASES),'jobs':list(jobs()),'source_hashes':source_hashes()}

def atomic_json(path,data):
    path=Path(path);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8');tmp.replace(path)
