import numpy as np
import scipy.linalg


'''for plotting'''
try:
    import pyvista
    pyvista.set_jupyter_backend('static') # Comment to enable interactive PyVista plots
    have_pyvista = True
except ImportError:
    have_pyvista = False


class VALENSystem:
    '''! Continuous-time valen system from Battey et al. equation (9)
    '''

    def __init__(self, s, a, tw_torus, tw_mode, coils=None):
        '''! Initialize a valen equation (9) state-space model

        @param Lw  wall-wall inductance matrix
        @param Lwc  wall-coil inductance matrix
        @param Lwd  wall-plasma inductance matrix
        @param Lcw  coil-wall inductance matrix
        @param Lc  coil-coil inductance matrix
        @param Lcd  coil-plasma inductance matrix
        @param Ldw  plasma-wall inductance matrix
        @param Ldc  plasma-coil inductance matrix
        @param Ld  plasma inductance matrix
        @param Rw Wall resistance matrix
        @param Rc coil-coil resistance matrix
        @param Rp Plasma-driver resistance matrix

        @param Lw_ef  Effective reluctance modified Lw
        @param Lwc_ef  Effective reluctance modified Lwc
        @param Lwd_ef  Effective reluctance modified Lwd
        @param Lcw_ef  Effective reluctance modified Lcw
        @param Lc_ef  Effective reluctance modified Lc
        @param Lcd_ef  Effective reluctance modified Lcd
        @param Ldw_ef  Effective reluctance modified Ldw
        @param Ldc_ef  Effective reluctance modified Ldc
        @param Ld_ef  Effective reluctance modified Ld
       
        @param a Boozer torque parameter 
        @param s Boozer stability parameter
        @param P Permiability matrix
        @param roe Reluctance matrix 
        '''

        self.tw_torus = tw_torus
        self.tw_mode = tw_mode
        
        self.tw_torus.compute_Lmat()
        '''how to compute without coils'''
        Lw = self.tw_torus.Lmat
        Mwd = self.tw_mode.cross_coupling(tw_torus)
        Mdw = Mwd
        self.tw_mode.compute_Lmat()
        '''Ld larger than Lw?'''
        Ld = self.tw_mode.Lmat

        self.tw_torus.compute_Rmat()
        Rw = self.tw_torus.Rmat
        self.tw_mode.compute_Rmat()
        Rd = self.tw_mode.Rmat
        

        if coils != None: 
            Mwc = self.tw_torus.compute_Mcoil()
            Mcw = Mwc
            Lc = self.tw_torus.n_icoils
            '''how to compute coil self inductance'''
            Mcd = self.tw_mode.compute_Mcoil()
            Mdc = Mcd 
            coils.compute_Rmat()
            Rc = coils.Rmat

        self.P = (-1/(s+a*1j))*np.identity(2)
        self.roe = np.linalg.inv(Ld)@(self.P-1)
        
        
        Lw_ef = Lw + Mwd @ self.roe @ Mdw
        Lwd_ef = Mwd + Mwd @ self.roe @ Ld
        Ldw_ef = Mdw + Ld @ self.roe @ Mdw
        Ld_ef = Ld + Ld @ self.roe@ Ld

        

        if coils != None: 
            Lwc_ef = Mwc + Mwd @ self.roe @ Mdc
            Ldc_ef = Mdc + Ld @ self.roe @ Mdc
            Lcw_ef = Mcw + Mcw @ self.roe @ Mdw 
            Lc_ef = Lc + Mcd @ self.roe @ Mcw
            Lcd_ef = Mcd + Mcd @ self.roe @ Ld 

            self.L = np.block([
                [Lw_ef, Lwc_ef, Lwd_ef],
                [Lcw_ef, Lc_ef, Lcd_ef],
                [Ldw_ef, Ldc_ef, Ld_ef]
            ])
            self.R = scipy.linalg.block_diag(Rw, Rc, Rd) 
        
        else: 
            self.L = np.block([
                [self.Lw, self.Lwd],
                [self.Ldw, self.Ld]
                ])
        
            self.R = scipy.linalg.block_diag(Rw, Rd)
        

    
    def rhs(self, currents, coil_voltage):
        '''! Evaluate \f$dI/dt\f$ for currents and coil voltages

        @param currents Full state vector \f$[I_w, I_c, I_d]^T\f$
        @param coil_voltage coil-coil voltage vector \f$V_c\f$
        @result Time derivative of the full state vector
        '''

        currents = np.asarray(currents, dtype=np.float64)
        coil_voltage = np.asarray(coil_voltage, dtype=np.float64)
        V = np.block([
            [np.zeros(len(coil_voltage))],
            [coil_voltage],
            [np.zeros(len(coil_voltage))]
            ])
        
        return scipy.linalg.solve(self.L, -np.dot(self.R, currents) + V)

    def eigenvalues(self):
        '''! Compute eigenvalues of the homogeneous equation (9) system

        @param with_vectors Return eigenvectors in addition to eigenvalues?
        '''
        #return scipy.linalg.eigs(self.L, -self.R), right=True)
        self.eig_vals, self.eig_vecs = scipy.sparse.linalg.eigs(self.L @ -self.R, k=2)
        return self.eig_vals, self.eig_vecs
    
    def plot(self):
        self.tw_mode.save_current(self.eig_vecs[0,:],'J_01')
        self.tw_mode.save_current(self.eig_vecs[1,:],'J_02')
        plot_data = self.tw_mode.build_XDMF()

        if have_pyvista:
            grid = plot_data['ThinCurr']['smesh'].get_pyvista_grid()
            J_01 = plot_data['ThinCurr']['smesh'].get_field('J_01_v')

            grid["vectors"] = J_01
            grid.set_active_vectors("vectors")

            p = pyvista.Plotter()
            scale = 0.1/(np.linalg.norm(J_01,axis=1)).max()
            arrows = grid.glyph(scale="vectors", orient="vectors", factor=scale)
            p.add_mesh(arrows, cmap="turbo", scalar_bar_args={'title': "Imag(J)", "vertical": True, "position_y":0.25, "position_x": 0.0})
            p.add_mesh(grid, color="white", opacity=1.0, show_edges=False)
            p.show(jupyter_backend='static')



    def plot(self):
        driver = np.zeros((2,self.tw_torus.nelems))
        driver[0,:] = Mcoil[0,:]*coil_current
        result = self.tw_torus.compute_freq_response(fdriver=driver,freq=1.E3)
        self.tw_torus.save_current(result[0,:],'Jr_coil')
        self.tw_torus.save_current(result[1,:],'Ji_coil')
        _ = self.tw_torus.build_XDMF()

        if have_pyvista:
            grid = plot_data['ThinCurr']['smesh'].get_pyvista_grid()
            Ji = plot_data['ThinCurr']['smesh'].get_field('Ji_mode_v')

            grid["vectors"] = Ji
            grid.set_active_vectors("vectors")

            p = pyvista.Plotter()
            scale = 0.2/(np.linalg.norm(Ji,axis=1)).max()
            arrows = grid.glyph(scale="vectors", orient="vectors", factor=scale)
            p.add_mesh(arrows, cmap="turbo", scalar_bar_args={'title': "Imag(J)", "vertical": True, "position_y":0.25, "position_x": 0.05})
            p.add_mesh(grid, color="white", opacity=1.0, show_edges=False)
            p.show(jupyter_backend='static')



    def flux(self):
        '''solve flux equation for certain eigen vals/vecs? '''
        
        Id = self.tw_mode.save_current(self.eig_vecs[0,:],'J_01')
        Ic = self.coil.save_current(self.eig_vecs[0,:],'J_01')
        '''coil model ??'''



# if no coils, remove coils from matrix 
# add in equations 2 and 3 
# solve for 'top' 2 eigenvalues 
# user only inputs model, class creates all ther variables 
# eigenvals frequency visualization? with currents? 

#Lc = coils.compute_Lmat() ???
#reduced thincurr model?? 


change i coil to v coil in xml file 
from git doc get first line of coil def with radius and ... 
use ex_oft_creation to make xml coil  
    put two circular coils outside torus? 
    could make c coils (two circular coils along the toroidal path) 
        circle with arbitrary x,y,z points 

use n_vcoils vairbale to get botton number of rows/colum values from block matrix, in bottom right ConnectionResetErrorLw = wall_model.Lmat[:-wall_model.nv_voils,:-wall_jodel.n_vcoils]

dont use compute_Mcoil(), thats for i coils 


plot eig vals on real-imaginary plane 
plot eignenvectors 

Lp and Pjj are both 2x2 



mode_driver = tw_mode.cross_eval(tw_wall,mode_driver) *** get code from an example 
mode_driver, 2x2, mutual for wall and plasma 

mode_driveT * self inductance for plasma * mode_driver
or maybe mode_drive * self inductance for plasma * mode_driverT

coil to plasma - dot current distrubution with big mutual matrix 

choose s and alpha from plots in battey (0 or smth for a, test diff s, should become unstable at zero) 