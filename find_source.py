import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from math import ceil

# ==============================================================================
# Dependency Check & Imports
# ==============================================================================
# Ensure pygacode can be imported
# Attempt to find pygacode in standard locations or relative to this script
sys.path.append(os.path.join(os.path.dirname(__file__), '../gacode/f2py'))

try:
    from pygacode.cgyro.data import cgyrodata
except ImportError:
    print("CRITICAL ERROR: Could not import 'pygacode.cgyro.data.cgyrodata'.")
    print("Please ensure pygacode is available in your PYTHONPATH or the script directory structure is correct.")
    sys.exit(1)

# ==============================================================================
# Class Definition: CGYRO_EnergyAnalysis
# ==============================================================================
class CGYRO_EnergyAnalysis(cgyrodata):
    """
    Extends cgyrodata with nonlinear energy transfer analysis capabilities.
    Adapted from OMFITcgyro_nonlin in collect.py.
    """
    def __init__(self, path, window=0.3, t_end=1.0):
        # Initialize the base cgyrodata class
        # cgyrodata expects the directory path
        super().__init__(path)
        print(f'Initialized CGYRO_EnergyAnalysis for {path}')
        
        # Ensure big field data is loaded
        print("Loading big field data (getbigfield)...")
        self.getbigfield()
        
        # Setup time averaging indices
        # self.t is available in cgyrodata
        if hasattr(self, 't') and len(self.t) > 0:
            self.n_time = len(self.t)
            t_max = self.t[-1]
            
            # Define window relative to t_end
            t_start_idx = int(self.n_time * (t_end - window))
            t_end_idx = int(self.n_time * t_end)
            
            if t_start_idx < 0: t_start_idx = 0
            if t_end_idx > self.n_time: t_end_idx = self.n_time
            
            self.ind_t_ave = np.arange(t_start_idx, t_end_idx)
            print(f"Time averaging window: indices {t_start_idx} to {t_end_idx} (t = {self.t[t_start_idx]:.2f} - {self.t[t_end_idx-1]:.2f})")
        else:
            print("Warning: No time data found.")
            self.ind_t_ave = np.array([])

        # Parse input.cgyro for IPCCW
        self.ipccw = 1.0
        try:
            input_file = os.path.join(path, 'input.cgyro')
            if os.path.exists(input_file):
                with open(input_file, 'r') as f:
                    for line in f:
                        if 'IPCCW' in line:
                            # Assumes format "IPCCW = 1" or similar
                            parts = line.split('=')
                            if len(parts) > 1:
                                val_str = parts[1].split()[0] # Take first token after =
                                self.ipccw = float(val_str)
                                print(f"Found IPCCW = {self.ipccw} in input.cgyro")
                                break
        except Exception as e:
            print(f"Warning: Could not parse IPCCW from input.cgyro: {e}")

        # Prepare complex fields for analysis
        if not hasattr(self, 'kxky_phi'):
            raise ValueError("kxky_phi not found in data.")

        shape = self.kxky_phi.shape
        ndim = self.kxky_phi.ndim
        
        if ndim == 5:
            # [2, n_radial, n_theta, n_n, n_time]
            self.phi_cmplx_t = self.kxky_phi[0, :, :, :, :] + 1j * self.kxky_phi[1, :, :, :, :]
        elif ndim == 6:
             print(f"Warning: kxky_phi has 6 dims: {shape}. Attempting to use index 0 for dim 3.")
             self.phi_cmplx_t = self.kxky_phi[0, :, :, 0, :, :] + 1j * self.kxky_phi[1, :, :, 0, :, :]
        elif ndim == 4:
            # [n_radial, n_theta, n_n, n_time] - likely already complex
            if np.iscomplexobj(self.kxky_phi):
                self.phi_cmplx_t = self.kxky_phi
            else:
                self.phi_cmplx_t = self.kxky_phi
        else:
            raise ValueError(f"Unexpected kxky_phi dimensions: {shape}")

    def Energy_transfer_phi(self, i_theta_plot=1, kx_select=0.0, ky_select=0.2):
        """
        Calculate kinetic energy transfer to a selected mode (kx_select, ky_select).
        """
        # Use actual data shape dimensions
        n_r_data = self.phi_cmplx_t.shape[0]
        n_n_data = self.phi_cmplx_t.shape[2]
        
        # Ensure i_theta_plot is within bounds
        if i_theta_plot >= self.phi_cmplx_t.shape[1]:
            i_theta_plot = 0
            print("Warning: i_theta_plot out of bounds, using 0.")

        # Normalize phi by rho to match get_phi_n behavior (phi_norm = phi_raw / rho)
        phi_cmplx = self.phi_cmplx_t[:, i_theta_plot, :, :] / self.rho
        
        # Find index of kx_select and ky_select
        kx_select_incode = self.kx.flat[np.abs(self.kx - kx_select).argmin()]
        ind_kx = np.where(self.kx == kx_select_incode)[0][0]
        
        ky_select_incode = self.ky.flat[np.abs(self.ky - ky_select).argmin()]
        ind_ky = np.where(self.ky == ky_select_incode)[0][0]
        
        # Calculation
        len_ind_t_ave = len(self.ind_t_ave)
        S_k_kp = np.zeros([n_r_data, n_n_data], dtype=complex)
        Lamda = np.zeros([n_r_data, n_n_data])
        ind_kx_0 = np.where(self.kx == 0)
        
        for i_m in np.arange(n_r_data): # Loop over kx'
            for i_n in np.arange(n_n_data): # Loop over ky'
                # Check if i_m or i_n are within bounds of kx/ky arrays
                if i_m >= len(self.kx) or i_n >= len(self.ky):
                    continue
                    
                kxp = self.kx[i_m]
                kyp = self.ky[i_n]
                kx_kxp = kx_select_incode - kxp
                ky_kyp = ky_select_incode - kyp
                
                # Coupling coefficient
                Lamda[i_m, i_n] = 1./2 * (kx_select_incode * kyp - kxp * ky_select_incode)
                Lamda[i_m, i_n] = Lamda[i_m, i_n] * ((kxp**2 + kyp**2) - (kx_kxp**2 + ky_kyp**2))
                
                S_k_kp_temp = np.zeros([len_ind_t_ave], dtype=complex)
                
                # Vectorized time loop for speed
                t_indices = self.ind_t_ave
                
                phi_k_select_t = phi_cmplx[ind_kx, ind_ky, t_indices]
                phi_k_p_t = phi_cmplx[i_m, i_n, t_indices]
                
                # Find phi(k-k')
                phi_k_kp_t = np.zeros(len_ind_t_ave, dtype=complex)
                
                if not (abs(ind_kx - i_m) > n_r_data // 2 - 1 or abs(ind_ky - i_n) > n_n_data):
                    ind_kx_m = ind_kx - i_m
                    ind_ky_kyp = 0 + (ind_ky - i_n)
                    
                    ind_kx_kxp = -1
                    use_conj = False
                    
                    if ind_ky_kyp < 0: # Conjugate symmetry
                        if len(ind_kx_0[0]) > 0:
                            ind_kx_kxp = ind_kx_0[0][0] - ind_kx_m
                            use_conj = True
                            target_ky_idx = -1 * ind_ky_kyp
                    else:
                        if len(ind_kx_0[0]) > 0:
                            ind_kx_kxp = ind_kx_0[0][0] + ind_kx_m
                            use_conj = False
                            target_ky_idx = ind_ky_kyp
                    
                    # Safe access
                    if ind_kx_kxp >= 0 and ind_kx_kxp < n_r_data and target_ky_idx >= 0 and target_ky_idx < n_n_data:
                        if use_conj:
                             phi_k_kp_t = np.conj(phi_cmplx[ind_kx_kxp, target_ky_idx, t_indices])
                        else:
                             phi_k_kp_t = phi_cmplx[ind_kx_kxp, target_ky_idx, t_indices]

                S_k_kp_temp = np.conj(phi_k_select_t) * phi_k_p_t * phi_k_kp_t
                S_k_kp[i_m, i_n] = np.real(np.mean(S_k_kp_temp))
        
        # Apply sign correction from IPCCW
        S_k_kp = -1 * self.ipccw * S_k_kp
        self.T_phi = Lamda * S_k_kp
        self.lamda = Lamda

    def get_transfer_spectrum(self, target_kx, target_ky, source_kx, source_ky, i_theta_plot=1):
        """
        Calculate the frequency spectrum of the energy transfer from source to target.
        """
        # Ensure i_theta_plot is within bounds
        if i_theta_plot >= self.phi_cmplx_t.shape[1]:
            i_theta_plot = 0
            # print("Warning: i_theta_plot out of bounds in get_transfer_spectrum, using 0.")

        # Find indices
        n_r_data = self.phi_cmplx_t.shape[0]
        n_n_data = self.phi_cmplx_t.shape[2]
        
        ind_kx = np.abs(self.kx - target_kx).argmin()
        ind_ky = np.abs(self.ky - target_ky).argmin()
        
        i_m = np.abs(self.kx - source_kx).argmin()
        i_n = np.abs(self.ky - source_ky).argmin()
        
        # Recalculate Lamda for this specific pair
        kx_select_incode = self.kx[ind_kx]
        ky_select_incode = self.ky[ind_ky]
        kxp = self.kx[i_m]
        kyp = self.ky[i_n]
        
        kx_kxp = kx_select_incode - kxp
        ky_kyp = ky_select_incode - kyp
        
        lamda = 1./2 * kx_select_incode * kyp - kxp * ky_select_incode
        lamda = lamda * ((kxp**2 + kyp**2) - (kx_kxp**2 + ky_kyp**2))
        
        # Time series
        t_indices = self.ind_t_ave
        # Normalize phi by rho to match get_phi_n behavior
        phi_cmplx = self.phi_cmplx_t[:, i_theta_plot, :, :] / self.rho
        
        phi_k_select_t = phi_cmplx[ind_kx, ind_ky, t_indices]
        phi_k_p_t = phi_cmplx[i_m, i_n, t_indices]
        
        # Find phi(k-k')
        phi_k_kp_t = np.zeros(len(t_indices), dtype=complex)
        ind_kx_0 = np.where(self.kx == 0)
        
        if not (abs(ind_kx - i_m) > n_r_data // 2 - 1 or abs(ind_ky - i_n) > n_n_data):
            ind_kx_m = ind_kx - i_m
            ind_ky_kyp = 0 + (ind_ky - i_n)
            ind_kx_kxp = -1
            use_conj = False
            target_ky_idx = -1
            
            if ind_ky_kyp < 0:
                if len(ind_kx_0[0]) > 0:
                    ind_kx_kxp = ind_kx_0[0][0] - ind_kx_m
                    use_conj = True
                    target_ky_idx = -1 * ind_ky_kyp
            else:
                if len(ind_kx_0[0]) > 0:
                    ind_kx_kxp = ind_kx_0[0][0] + ind_kx_m
                    use_conj = False
                    target_ky_idx = ind_ky_kyp
            
            if ind_kx_kxp >= 0 and ind_kx_kxp < n_r_data and target_ky_idx >= 0 and target_ky_idx < n_n_data:
                if use_conj:
                    phi_k_kp_t = np.conj(phi_cmplx[ind_kx_kxp, target_ky_idx, t_indices])
                else:
                    phi_k_kp_t = phi_cmplx[ind_kx_kxp, target_ky_idx, t_indices]

        # Calculate Transfer Time Series
        # S(t) = conj(phi_k)*phi_p*phi_q
        # Calculate Transfer Time Series
        # S(t) = conj(phi_k)*phi_p*phi_q
        # phi is already normalized by rho above
        S_t = np.conj(phi_k_select_t) * phi_k_p_t * phi_k_kp_t
        
        T_phi_t = -1 * self.ipccw * lamda * np.real(S_t) # Take real part for energy transfer?
        # Usually Energy transfer is Re(S_t).
        
        # Perform FFT
        dt = self.t[1] - self.t[0]
        omega = np.fft.fftfreq(len(T_phi_t), d=dt) * 2 * np.pi
        T_omega = np.fft.fft(T_phi_t)
        
        # Shift
        omega_shifted = np.fft.fftshift(omega)
        T_omega_shifted = np.fft.fftshift(T_omega)
        
        return omega_shifted, np.abs(T_omega_shifted)

# ==============================================================================
# Tracking Algorithm
# ==============================================================================
def track_energy_source(case_obj, start_kx, start_ky, max_steps=20, theta_idx=None):
    """
    Trace the energy source backwards from (start_kx, start_ky).
    """
    if theta_idx is None:
        theta_idx = case_obj.theta_plot // 2

    path = []
    current_kx = start_kx
    current_ky = start_ky

    print(f"  [Start Tracking] Target: ({current_kx:.3f}, {current_ky:.3f})")
    
    for step in range(max_steps):
        path.append((current_kx, current_ky))
        
        # 1. Calculate Energy Transfer to Current Target
        try:
            case_obj.Energy_transfer_phi(i_theta_plot=theta_idx, kx_select=current_kx, ky_select=current_ky)
        except Exception as e:
            print(f"    Error at step {step}: {e}")
            break
            
        # 2. Find Source Mode with Max Positive Transfer
        # T_phi > 0 means Energy flows FROM (kx', ky') TO (current_kx, current_ky)
        T_matrix = case_obj.T_phi
        lamda = case_obj.lamda
        
        # Find max index
        max_idx_flat = np.argmax(T_matrix)
        idx_kx_src, idx_ky_src = np.unravel_index(max_idx_flat, T_matrix.shape)
        
        src_kx = case_obj.kx[idx_kx_src]
        src_ky = case_obj.ky[idx_ky_src]
        max_val = T_matrix[idx_kx_src, idx_ky_src]
        
        print(f"    Step {step+1}: Target({current_kx:.3f}, {current_ky:.3f}) <== Source({src_kx:.3f}, {src_ky:.3f}) [T={max_val:.2e}] [Lambda={lamda[idx_kx_src, idx_ky_src]}]")
        
        # 3. Checks
        # Stop if max transfer is negative or zero (no driving source)
        if max_val <= 0:
            print(f"    --> STOP: Max energy transfer is non-positive ({max_val:.2e}). No driving source found.")
            # path.append((src_kx, src_ky)) # Do not append, as this source is not valid
            break

        # Loop detection
        if (src_kx, src_ky) in path:
            print(f"    --> Loop detected! ({src_kx:.3f}, {src_ky:.3f}) visited before.")
            path.append((src_kx, src_ky)) # Add closing point
            break
            
        # Convergence to kx ~ 0
        if abs(src_kx) < 1e-4:
            print(f"    --> Converged to kx~0 at ky={src_ky:.3f}")
            path.append((src_kx, src_ky))
            break
            
        # Update for next step
        current_kx = src_kx
        current_ky = src_ky
    
    return path

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    # Check for arguments
    if len(sys.argv) > 1:
        sim_dir = sys.argv[1]
    else:
        sim_dir = '.' # Current directory
        
    print(f"Analyzing simulation in: {os.path.abspath(sim_dir)}")
    
    try:
        case = CGYRO_EnergyAnalysis(sim_dir)
    except Exception as e:
        print(f"Failed to load CGYRO data: {e}")
        import traceback
        traceback.print_exc()
        return

    # Algorithm Settings
    available_kys = case.ky
    # Filter out ky=0 or extremely small ky
    start_kys = [ky for ky in available_kys if ky > 1e-4]
    
    print(f"Scanning {len(start_kys)} ky modes: {start_kys}")

    # Plotting Setup
    # Figure 1: Paths
    fig1 = plt.figure(figsize=(10, 8))
    ax1 = fig1.add_subplot(111)
    
    kx_grid, ky_grid = np.meshgrid(case.kx, case.ky)
    ax1.scatter(kx_grid.flatten(), ky_grid.flatten(), c='lightgray', s=10, alpha=0.5, label='Grid')

    # Figure 2: FFT Spectra (Subplots)
    # We will create a separate figure for spectra, or save individual plots?
    # Given potentially many starting kys, plotting all on one might be messy.
    # Let's create one figure with many subplots or just plot the most significant ones.
    # For now, let's create a separate PDF or just save individual PNGs if requested.
    # Or, plot them all on one plot with offset?
    
    # Let's create a second figure for the spectra
    fig2 = plt.figure(figsize=(12, 6))
    ax2 = fig2.add_subplot(111)
    ax2.set_title("FFT of Energy Transfer for Converged Paths")
    ax2.set_xlabel(r'$\omega (c_s/a)$')
    ax2.set_ylabel(r'$|FFT(T_\phi)|$')

    # Use a colormap that can handle many lines
    cmap = plt.cm.jet
    
    for i, ky_start in enumerate(start_kys):
        real_ky_start = ky_start
        real_kx_start = 0.0 # Start from kx=0
        
        print(f"\n--- Tracking for Start ky = {real_ky_start:.3f} ---")
        
        path = track_energy_source(case, real_kx_start, real_ky_start)
        
        # Visualize Path
        path = np.array(path)
        color = cmap(float(i) / len(start_kys))
        
        if len(path) > 1:
            # Draw line (thinner)
            ax1.plot(path[:, 0], path[:, 1], '-', color=color, linewidth=1.0, alpha=0.8)
            
            # Draw intermediate points (small)
            if len(path) > 2:
                ax1.plot(path[1:-1, 0], path[1:-1, 1], '.', color=color, markersize=4)
            
            # Mark Start (larger)
            ax1.plot(path[0, 0], path[0, 1], 'o', color=color, markersize=6, label=f'ky={real_ky_start:.2f}')
            # Mark End (larger)
            ax1.plot(path[-1, 0], path[-1, 1], 's', color=color, markersize=8, markeredgewidth=1.5, markerfacecolor='none')
            
            # --- FFT Analysis of the Final Step ---
            # We want the transfer INTO the last mode in the path FROM the mode before it?
            # Or the transfer INTO the second-to-last mode FROM the last mode?
            # The tracking finds Source -> Target.
            # path[-1] is the Source found in the last step.
            # path[-2] is the Target in the last step.
            # So we calculate T(Target <- Source)
            
            target_kx, target_ky = path[-2]
            source_kx, source_ky = path[-1]
            
            print(f"    Calculating FFT for Transfer: ({source_kx:.2f},{source_ky:.2f}) -> ({target_kx:.2f},{target_ky:.2f})")
            
            omega, t_spec = case.get_transfer_spectrum(target_kx, target_ky, source_kx, source_ky)
            
            # Normalize spectrum for comparison
            if np.max(t_spec) > 0:
                t_spec = t_spec / np.max(t_spec)
                
            ax2.plot(omega, t_spec, color=color, alpha=0.6, linewidth=1.0) #, label=f'ky={real_ky_start:.2f}')


    ax1.set_xlabel(r'$k_x \rho_s$')
    ax1.set_ylabel(r'$k_y \rho_s$')
    ax1.set_title('Reverse Energy Transfer Tracking')
    ax1.grid(True, alpha=0.3)
    # ax1.legend() # Legend might be too big if many modes
    
    output_plot_path = 'energy_tracking_path.png'
    #fig1.savefig(output_plot_path, dpi=150)
    print(f"Path plot saved to {output_plot_path}")
    
    output_plot_fft = 'energy_transfer_fft.png'
    #fig2.savefig(output_plot_fft, dpi=150)
    print(f"FFT plot saved to {output_plot_fft}")
    
    plt.show()

if __name__ == "__main__":
    main()

