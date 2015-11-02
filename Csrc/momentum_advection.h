#pragma once
#include "grid.h"
#include "advection_interpolation.h"
#include "flux_divergence.h"
#include<stdio.h>

void second_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 0;
        const ssize_t jmin = 0;
        const ssize_t kmin = 0;

        const ssize_t imax = dims->nlg[0]-1;
        const ssize_t jmax = dims->nlg[1]-1;
        const ssize_t kmax = dims->nlg[2]-1;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp1_ing = stencil[d_advected];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) )*rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) )*rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) )*rho0[k];
                    }
                }
            }
        }

        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);

        //Free dynamically allocated array
        free(flux);
        return;
    }


void second_order_ml_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 0;
        const ssize_t jmin = 0;
        const ssize_t kmin = 0;

        const ssize_t imax = dims->nlg[0]-1;
        const ssize_t jmax = dims->nlg[1]-1;
        const ssize_t kmax = dims->nlg[2]-1;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp1_ing = stencil[d_advected];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) )*rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk]*rho0[k],vel_advecting[ijk+sp1_ing]*rho0[k])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) );
                    }
                }
            }
        }
        else if(d_advected == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk] * rho0_half[k],vel_advecting[ijk+sp1_ing]*rho0_half[k+1])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) );
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing])
                            *interp_2(vel_advected[ijk],vel_advected[ijk + sp1_ed]) )*rho0[k];
                    }
                }
            }
        }

        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);

        //Free dynamically allocated array
        free(flux);
        return;
    }

void fourth_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
        double* restrict alpha0, double* restrict alpha0_half,
        double* restrict vel_advected, double* restrict vel_advecting,
        double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);


        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 1;
        const ssize_t jmin = 1;
        const ssize_t kmin = 1;

        const ssize_t imax = dims->nlg[0]-2;
        const ssize_t jmax = dims->nlg[1]-2;
        const ssize_t kmax = dims->nlg[2]-2;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_4(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_4(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_4(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed]))* rho0[k];
                    }
                }
            }
        }

        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);

        free(flux);
        return;
    }

void fourth_order_ws_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 1;
        const ssize_t jmin = 1;
        const ssize_t kmin = 1;

        const ssize_t imax = dims->nlg[0]-2;
        const ssize_t jmax = dims->nlg[1]-2;
        const ssize_t kmax = dims->nlg[2]-2;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_4(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0[k];
                    }
                }
            }
        }

        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void fourth_order_ml_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Y. Morinishi, T.S. Lund, O.V. Vasilyev, P. Moin, Fully Conservative Higher Order Finite Difference Schemes
        // for Incompressible Flow, Journal of Computational Physics, Volume 143, Issue 1, 10 June 1998,
        // Pages 90-124, ISSN 0021-9991, http://dx.doi.org/10.1006/jcph.1998.5962.

        // Dynamically allocate flux array
        double *flux1 = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);
        double *flux3 = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);


        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 1;
        const ssize_t jmin = 1;
        const ssize_t kmin = 1;

        const ssize_t imax = dims->nlg[0]-2;
        const ssize_t jmax = dims->nlg[1]-2;
        const ssize_t kmax = dims->nlg[2]-2;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = 2 * sm1_ed;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;

        const double dxi = dims->dxi[d_advecting];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){

                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (9.0/16.0*( vel_advecting[ijk] + vel_advecting[ijk + sp1_ing])
                                          - 1.0/16.0*(vel_advecting[ijk+sm1_ing] + vel_advecting[ijk + sp2_ing]))*rho0_half[k];

                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])/2.0 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])/2.0 ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (9.0/16.0*( vel_advecting[ijk] * rho0[k] +  vel_advecting[ijk + sp1_ing]*rho0[k+1])
                                           - 1.0/16.0*(vel_advecting[ijk+sm1_ing]*rho0[k-1] + vel_advecting[ijk + sp2_ing]*rho0[k+2]));
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])/2.0 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])/2.0 ;
                    }
                }
            }
        }
        else if(d_advected == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (9.0/16.0*( vel_advecting[ijk]*rho0_half[k] +  vel_advecting[ijk + sp1_ing]*rho0_half[k+1])
                                           - 1.0/16.0*(vel_advecting[ijk+sm1_ing]*rho0_half[k-1] + vel_advecting[ijk + sp2_ing]*rho0_half[k+2] ));
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])/2.0 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])/2.0 ;
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (9.0/16.0*(vel_advecting[ijk] +  vel_advecting[ijk + sp1_ing])
                                           - 1.0/16.0*(vel_advecting[ijk+sm1_ing] + vel_advecting[ijk + sp2_ing]))*rho0[k];
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])/2.0 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])/2.0 ;
                    }
                }
            }
        }
        /// Now compute the velocity tendencies
        if(d_advected != 2){
            for(ssize_t i=imin + 1; i<imax; i++){
                const ssize_t ishift = i * istride;
                for(ssize_t j=jmin + 1; j<jmax; j++){
                    const ssize_t jshift = j * jstride;
                    for(ssize_t k=kmin + 1 ; k<kmax; k++){
                        const ssize_t ijk = ishift + jshift + k;
                        tendency[ijk] -= alpha0_half[k] * (9.0 /16.0 * (flux1[ijk] - flux1[ijk + sm1_ed])
                                                             - 1.0/24.0 * (flux3[ijk+sp1_ed] - flux3[ijk + sm2_ed]) )*dxi;
                    }
                 }
            }
        }
        else{
            for(ssize_t i=imin + 1 ; i<imax; i++){
                const ssize_t ishift = i * istride;
                for(ssize_t j=jmin + 1 ; j<jmax; j++){
                    const ssize_t jshift = j * jstride;
                    for(ssize_t k=kmin + 1 ; k<kmax; k++){
                        const ssize_t ijk = ishift + jshift + k;
                        tendency[ijk] -= alpha0[k] * (9.0 /16.0 * (flux1[ijk] - flux1[ijk + sm1_ed])
                                                             - 1.0/24.0 * (flux3[ijk+sp1_ed] - flux3[ijk + sm2_ed]) )*dxi;
                    }
                 }
            }
        }

        free(flux1);
        free(flux3);
        return;
    }



void sixth_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){


        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 2;
        const ssize_t jmin = 2;
        const ssize_t kmin = 2;

        const ssize_t imax = dims->nlg[0]-3;
        const ssize_t jmax = dims->nlg[1]-3;
        const ssize_t kmax = dims->nlg[2]-3;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sp3_ed = 3 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = -2*sp1_ed;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sp3_ing = 3 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;
        const ssize_t sm2_ing = -2*sp1_ing;

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_6(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed])) * rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_6(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed]))* rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_6(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed]))* rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void sixth_order_ws_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 2;
        const ssize_t jmin = 2;
        const ssize_t kmin = 2;

        const ssize_t imax = dims->nlg[0]-3;
        const ssize_t jmax = dims->nlg[1]-3;
        const ssize_t kmax = dims->nlg[2]-3;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sp3_ed = 3 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = -2*sp1_ed;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sp3_ing = 3 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;
        const ssize_t sm2_ing = -2*sp1_ing;

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed])) * rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed])) * rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_2(vel_advecting[ijk],vel_advecting[ijk+sp1_ing]) *
                                 interp_6(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed])) * rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void sixth_order_ml_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Y. Morinishi, T.S. Lund, O.V. Vasilyev, P. Moin, Fully Conservative Higher Order Finite Difference Schemes
        // for Incompressible Flow, Journal of Computational Physics, Volume 143, Issue 1, 10 June 1998,
        // Pages 90-124, ISSN 0021-9991, http://dx.doi.org/10.1006/jcph.1998.5962.

        // Dynamically allocate flux array
        double *flux1 = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);
        double *flux3 = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);
        double *flux5 = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 2;
        const ssize_t jmin = 2;
        const ssize_t kmin = 2;

        const ssize_t imax = dims->nlg[0]-3;
        const ssize_t jmax = dims->nlg[1]-3;
        const ssize_t kmax = dims->nlg[2]-3;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sp3_ed = 3 * sp1_ed;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = 2 * sm1_ed;
        const ssize_t sm3_ed = 3 * sm1_ed;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sp3_ing = 3 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;
        const ssize_t sm2_ing = -2 * sp1_ing;


        const double dxi = dims->dxi[d_advecting];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){

                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (150.0/256.0*( vel_advecting[ijk] + vel_advecting[ijk + sp1_ing])
                                          - 25.0/256.0*(vel_advecting[ijk+sm1_ing] + vel_advecting[ijk + sp2_ing])
                                          + 3.0/256.0*(vel_advecting[ijk+sm2_ing] + vel_advecting[ijk + sp3_ing]))*rho0_half[k];

                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])*0.5 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])*0.5 ;
                        flux5[ijk] = d * (vel_advected[ijk+sm2_ed] + vel_advected[ijk+sp3_ed])*0.5 ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (150.0/256.0*( vel_advecting[ijk]*rho0[k] + vel_advecting[ijk + sp1_ing]*rho0[k+1])
                                          - 25.0/256.0*(vel_advecting[ijk+sm1_ing]*rho0[k-1] + vel_advecting[ijk + sp2_ing]*rho0[k+2])
                                          + 3.0/256.0*(vel_advecting[ijk+sm2_ing]*rho0[k-2] + vel_advecting[ijk + sp3_ing]*rho0[k+3]));
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])*0.5 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])*0.5 ;
                        flux5[ijk] = d * (vel_advected[ijk+sm2_ed] + vel_advected[ijk+sp3_ed])*0.5 ;
                    }
                }
            }
        }
        else if(d_advected == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (150.0/256.0*( vel_advecting[ijk]*rho0_half[k] + vel_advecting[ijk + sp1_ing]*rho0_half[k+1])
                                          - 25.0/256.0*(vel_advecting[ijk+sm1_ing]*rho0_half[k-1] + vel_advecting[ijk + sp2_ing]*rho0_half[k+2])
                                          + 3.0/256.0*(vel_advecting[ijk+sm2_ing]*rho0_half[k-2] + vel_advecting[ijk + sp3_ing]*rho0_half[k+3]));
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])*0.5 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])*0.5 ;
                        flux5[ijk] = d * (vel_advected[ijk+sm2_ed] + vel_advected[ijk+sp3_ed])*0.5 ;
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        const double d = (150.0/256.0*(rho0[k]*vel_advecting[ijk] +  rho0[k+1]*vel_advecting[ijk + sp1_ing])
                                           - 25.0/256.0*(rho0[k-1]*vel_advecting[ijk+sm1_ing] + rho0[k+2]*vel_advecting[ijk + sp2_ing])
                                            + 3.0/256.0*(rho0[k-2]*vel_advecting[ijk+sm2_ing] + rho0[k+3]*vel_advecting[ijk + sp3_ing]))*rho0[k];
                        flux1[ijk] = d * (vel_advected[ijk] +  vel_advected[ijk + sp1_ed])*0.5 ;
                        flux3[ijk] = d * (vel_advected[ijk+sm1_ed] + vel_advected[ijk+sp2_ed])*0.5 ;
                        flux5[ijk] = d * (vel_advected[ijk+sm2_ed] + vel_advected[ijk+sp3_ed])*0.5 ;
                    }
                }
            }
        }

        /// Now compute the velocity tendencies

        if(d_advected != 2){
            for(ssize_t i=imin + 3; i<imax-2; i++){
                const ssize_t ishift = i * istride;
                for(ssize_t j=jmin + 3; j<jmax-2; j++){
                    const ssize_t jshift = j * jstride;
                    for(ssize_t k=kmin + 3 ; k<kmax-2; k++){
                        const ssize_t ijk = ishift + jshift + k;
                        tendency[ijk] -= alpha0_half[k] * (150.0/256.0 * (flux1[ijk] - flux1[ijk + sm1_ed])
                                                             - 25.0/384.0 * (flux3[ijk+sp1_ed] - flux3[ijk + sm2_ed])
                                                              + 3.0/640.0 * (flux5[ijk+sp2_ed] - flux5[ijk + sm3_ed]))*dxi;
                    }
                 }
            }
        }
        else{
            for(ssize_t i=imin + 3 ; i<imax-2; i++){
                const ssize_t ishift = i * istride;
                for(ssize_t j=jmin + 3 ; j<jmax-2; j++){
                    const ssize_t jshift = j * jstride;
                    for(ssize_t k=kmin + 3 ; k<kmax-2; k++){
                        const ssize_t ijk = ishift + jshift + k;
                        tendency[ijk] -= alpha0[k] * (150.0 /256.0 * (flux1[ijk] - flux1[ijk + sm1_ed])
                                                             - 25.0/384.0 * (flux3[ijk+sp1_ed] - flux3[ijk + sm2_ed])
                                                              + 3.0/640.0 * (flux5[ijk+sp2_ed] - flux5[ijk + sm3_ed]))*dxi;
                    }
                 }
            }
        }

        free(flux1);
        free(flux3);
        free(flux5);
        return;
    }


void eighth_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 3;
        const ssize_t jmin = 3;
        const ssize_t kmin = 3;

        const ssize_t imax = dims->nlg[0]-4;
        const ssize_t jmax = dims->nlg[1]-4;
        const ssize_t kmax = dims->nlg[2]-4;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2 * sp1_ed ;
        const ssize_t sp3_ed = 3 * sp1_ed ;
        const ssize_t sp4_ed = 4 * sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = -2*sp1_ed;
        const ssize_t sm3_ed = -3*sp1_ed;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2 * sp1_ing;
        const ssize_t sp3_ing = 3 * sp1_ing;
        const ssize_t sp4_ing = 4 * sp1_ing;
        const ssize_t sm1_ing = -sp1_ing;
        const ssize_t sm2_ing = -2*sp1_ing;
        const ssize_t sm3_ing = -3*sp1_ing;

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_8(vel_advecting[ijk+sm3_ing],
                                          vel_advecting[ijk+sm2_ing],
                                          vel_advecting[ijk+sm1_ing],
                                          vel_advecting[ijk],
                                          vel_advecting[ijk+sp1_ing],
                                          vel_advecting[ijk+sp2_ing],
                                          vel_advecting[ijk+sp3_ing],
                                          vel_advecting[ijk+sp4_ing])
                                 *interp_8(vel_advected[ijk+sm3_ed],
                                           vel_advected[ijk+sm2_ed],
                                           vel_advected[ijk+sm1_ed],
                                           vel_advected[ijk],
                                           vel_advected[ijk+sp1_ed],
                                           vel_advected[ijk+sp2_ed],
                                           vel_advected[ijk+sp3_ed],
                                           vel_advected[ijk+sp4_ed])) * rho0_half[k];
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_8(vel_advecting[ijk+sm3_ing],
                                          vel_advecting[ijk+sm2_ing],
                                          vel_advecting[ijk+sm1_ing],
                                          vel_advecting[ijk],
                                          vel_advecting[ijk+sp1_ing],
                                          vel_advecting[ijk+sp2_ing],
                                          vel_advecting[ijk+sp3_ing],
                                          vel_advecting[ijk+sp4_ing])
                                 *interp_8(vel_advected[ijk+sm3_ed],
                                           vel_advected[ijk+sm2_ed],
                                           vel_advected[ijk+sm1_ed],
                                           vel_advected[ijk],
                                           vel_advected[ijk+sp1_ed],
                                           vel_advected[ijk+sp2_ed],
                                           vel_advected[ijk+sp3_ed],
                                           vel_advected[ijk+sp4_ed])) * rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                    flux[ijk] = (interp_8(vel_advecting[ijk+sm3_ing],
                                          vel_advecting[ijk+sm2_ing],
                                          vel_advecting[ijk+sm1_ing],
                                          vel_advecting[ijk],
                                          vel_advecting[ijk+sp1_ing],
                                          vel_advecting[ijk+sp2_ing],
                                          vel_advecting[ijk+sp3_ing],
                                          vel_advecting[ijk+sp4_ing])
                                 *interp_8(vel_advected[ijk+sm3_ed],
                                           vel_advected[ijk+sm2_ed],
                                           vel_advected[ijk+sm1_ed],
                                           vel_advected[ijk],
                                           vel_advected[ijk+sp1_ed],
                                           vel_advected[ijk+sp2_ed],
                                           vel_advected[ijk+sp3_ed],
                                           vel_advected[ijk+sp4_ed])) * rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void weno_third_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 1;
        const ssize_t jmin = 1;
        const ssize_t kmin = 1;

        const ssize_t imax = dims->nlg[0]-2;
        const ssize_t jmax = dims->nlg[1]-2;
        const ssize_t kmax = dims->nlg[2]-2;

        const ssize_t stencil[3] = {istride,jstride,1};
        const ssize_t sm1_ed = -stencil[d_advecting];
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2*stencil[d_advecting];
        const ssize_t sp1_ing = stencil[d_advected];


        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        //Upwind for positive velocity
                        const double phip = interp_weno3(vel_advected[ijk + sm1_ed],
                                                        vel_advected[ijk],
                                                        vel_advected[ijk + sp1_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno3(vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk]);

                        const double vel_adv = interp_2(vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        //Upwind for positive velocity
                        const double phip = interp_weno3(vel_advected[ijk + sm1_ed],
                                            vel_advected[ijk],
                                            vel_advected[ijk + sp1_ed]);

                        // Up wind for negative velocity
                        const double phim = interp_weno3(vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk]);
                        const double vel_adv = interp_2(vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        //Upwind for positive velocity
                        const double phip = interp_weno3(vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed]);

                        // Up wind for negative velocity
                        const double phim = interp_weno3(vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk]);

                        const double vel_adv = interp_2(vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void weno_fifth_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 2;
        const ssize_t jmin = 2;
        const ssize_t kmin = 2;

        const ssize_t imax = dims->nlg[0]-3;
        const ssize_t jmax = dims->nlg[1]-3;
        const ssize_t kmax = dims->nlg[2]-3;

        const ssize_t stencil[3] = {istride,jstride,1};

        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2*sp1_ed ;
        const ssize_t sp3_ed = 3*sp1_ed ;
        const ssize_t sm1_ed = -sp1_ed ;
        const ssize_t sm2_ed = -2*sp1_ed ;

        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2*sp1_ing ;
        const ssize_t sm1_ing = -sp1_ing ;


        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed]);

                        const double vel_adv = interp_4(vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed]);

                        const double vel_adv = interp_4(vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed]);

                        const double vel_adv = interp_4(vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void weno_seventh_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 3;
        const ssize_t jmin = 3;
        const ssize_t kmin = 3;

        const ssize_t imax = dims->nlg[0]-4;
        const ssize_t jmax = dims->nlg[1]-4;
        const ssize_t kmax = dims->nlg[2]-4;

        const ssize_t stencil[3] = {istride,jstride,1};

        const ssize_t sm3_ed = -3*stencil[d_advecting];
        const ssize_t sm2_ed = -2*stencil[d_advecting];
        const ssize_t sm1_ed = -stencil[d_advecting];
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2*stencil[d_advecting];
        const ssize_t sp3_ed = 3*stencil[d_advecting];
        const ssize_t sp4_ed = 4*stencil[d_advecting];

        const ssize_t sm2_ing = -2*stencil[d_advected];
        const ssize_t sm1_ing = -stencil[d_advected];
        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2*stencil[d_advected];
        const ssize_t sp3_ing = 3*stencil[d_advected];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp3_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed]);

                        const double vel_adv = interp_6(vel_advecting[ijk + sm2_ing],
                                                        vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp3_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed]);

                        const double vel_adv = interp_6(vel_advecting[ijk + sm2_ing],
                                                        vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;
                        //Upwind for positive velocity
                        const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp3_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed]);

                        const double vel_adv = interp_6(vel_advecting[ijk + sm2_ing],
                                                        vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }


void weno_ninth_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 4;
        const ssize_t jmin = 4;
        const ssize_t kmin = 4;

        const ssize_t imax = dims->nlg[0]-5;
        const ssize_t jmax = dims->nlg[1]-5;
        const ssize_t kmax = dims->nlg[2]-5;

        const ssize_t stencil[3] = {istride,jstride,1};

        const ssize_t sm4_ed = -4*stencil[d_advecting];
        const ssize_t sm3_ed = -3*stencil[d_advecting];
        const ssize_t sm2_ed = -2*stencil[d_advecting];
        const ssize_t sm1_ed = -stencil[d_advecting];
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2*stencil[d_advecting];
        const ssize_t sp3_ed = 3*stencil[d_advecting];
        const ssize_t sp4_ed = 4*stencil[d_advecting];
        const ssize_t sp5_ed = 5*stencil[d_advecting];

        const ssize_t sm3_ing = -3*stencil[d_advected];
        const ssize_t sm2_ing = -2*stencil[d_advected];
        const ssize_t sm1_ing = -stencil[d_advected];
        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2*stencil[d_advected];
        const ssize_t sp3_ing = 3*stencil[d_advected];
        const ssize_t sp4_ing = 4*stencil[d_advected];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                        vel_advected[ijk + sm3_ed],
                                                        vel_advected[ijk + sm2_ed],
                                                        vel_advected[ijk + sm1_ed],
                                                        vel_advected[ijk],
                                                        vel_advected[ijk + sp1_ed],
                                                        vel_advected[ijk + sp2_ed],
                                                        vel_advected[ijk + sp3_ed],
                                                        vel_advected[ijk + sp4_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                         vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm3_ed]);

                        const double vel_adv = interp_8(vel_advecting[ijk + sm3_ing],
                                                        vel_advecting[ijk + sm2_ing],
                                                        vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing],
                                                        vel_advecting[ijk + sp4_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;


                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                        vel_advected[ijk + sm3_ed],
                                                        vel_advected[ijk + sm2_ed],
                                                        vel_advected[ijk + sm1_ed],
                                                        vel_advected[ijk],
                                                        vel_advected[ijk + sp1_ed],
                                                        vel_advected[ijk + sp2_ed],
                                                        vel_advected[ijk + sp3_ed],
                                                        vel_advected[ijk + sp4_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                         vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm3_ed]);

                        const double vel_adv = interp_8(vel_advecting[ijk+sm3_ing],
                                                        vel_advecting[ijk+sm2_ing],
                                                        vel_advecting[ijk+sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing],
                                                        vel_advecting[ijk + sp4_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];

                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                        vel_advected[ijk + sm3_ed],
                                                        vel_advected[ijk + sm2_ed],
                                                        vel_advected[ijk + sm1_ed],
                                                        vel_advected[ijk],
                                                        vel_advected[ijk + sp1_ed],
                                                        vel_advected[ijk + sp2_ed],
                                                        vel_advected[ijk + sp3_ed],
                                                        vel_advected[ijk + sp4_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                         vel_advected[ijk + sp4_ed],
                                                         vel_advected[ijk + sp3_ed],
                                                         vel_advected[ijk + sp2_ed],
                                                         vel_advected[ijk + sp1_ed],
                                                         vel_advected[ijk],
                                                         vel_advected[ijk + sm1_ed],
                                                         vel_advected[ijk + sm2_ed],
                                                         vel_advected[ijk + sm3_ed]);

                        const double vel_adv = interp_8(vel_advecting[ijk + sm3_ing],
                                                        vel_advecting[ijk + sm2_ing],
                                                        vel_advecting[ijk + sm1_ing],
                                                        vel_advecting[ijk],
                                                        vel_advecting[ijk + sp1_ing],
                                                        vel_advecting[ijk + sp2_ing],
                                                        vel_advecting[ijk + sp3_ing],
                                                        vel_advecting[ijk + sp4_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];

                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

void weno_eleventh_order_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
    double* restrict alpha0, double* restrict alpha0_half,
    double* restrict vel_advected, double* restrict vel_advecting,
    double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

        const ssize_t istride = dims->nlg[1] * dims->nlg[2];
        const ssize_t jstride = dims->nlg[2];

        const ssize_t imin = 5;
        const ssize_t jmin = 5;
        const ssize_t kmin = 5;

        const ssize_t imax = dims->nlg[0]-6;
        const ssize_t jmax = dims->nlg[1]-6;
        const ssize_t kmax = dims->nlg[2]-6;

        const ssize_t stencil[3] = {istride,jstride,1};

        const ssize_t sm5_ed = -5*stencil[d_advecting];
        const ssize_t sm4_ed = -4*stencil[d_advecting];
        const ssize_t sm3_ed = -3*stencil[d_advecting];
        const ssize_t sm2_ed = -2*stencil[d_advecting];
        const ssize_t sm1_ed = -stencil[d_advecting];
        const ssize_t sp1_ed = stencil[d_advecting];
        const ssize_t sp2_ed = 2*stencil[d_advecting];
        const ssize_t sp3_ed = 3*stencil[d_advecting];
        const ssize_t sp4_ed = 4*stencil[d_advecting];
        const ssize_t sp5_ed = 5*stencil[d_advecting];
        const ssize_t sp6_ed = 6*stencil[d_advecting];

        const ssize_t sm4_ing = -4*stencil[d_advected];
        const ssize_t sm3_ing = -3*stencil[d_advected];
        const ssize_t sm2_ing = -2*stencil[d_advected];
        const ssize_t sm1_ing = -stencil[d_advected];
        const ssize_t sp1_ing = stencil[d_advected];
        const ssize_t sp2_ing = 2*stencil[d_advected];
        const ssize_t sp3_ing = 3*stencil[d_advected];
        const ssize_t sp4_ing = 4*stencil[d_advected];
        const ssize_t sp5_ing = 5*stencil[d_advected];

        if (d_advected != 2 && d_advecting !=2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno11(vel_advected[ijk + sm5_ed],
                                                          vel_advected[ijk + sm4_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp5_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno11(vel_advected[ijk + sp6_ed],
                                                          vel_advected[ijk + sp5_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm4_ed]);

                        const double vel_adv = interp_10(vel_advecting[ijk + sm4_ing],
                                                         vel_advecting[ijk + sm3_ing],
                                                         vel_advecting[ijk + sm2_ing],
                                                         vel_advecting[ijk + sm1_ing],
                                                         vel_advecting[ijk],
                                                         vel_advecting[ijk + sp1_ing],
                                                         vel_advecting[ijk + sp2_ing],
                                                         vel_advecting[ijk + sp3_ing],
                                                         vel_advecting[ijk + sp4_ing],
                                                         vel_advecting[ijk + sp5_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                    }
                }
            }
        }
        else if(d_advected == 2 && d_advecting == 2){
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno11(vel_advected[ijk + sm5_ed],
                                                          vel_advected[ijk + sm4_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp5_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno11(vel_advected[ijk + sp6_ed],
                                                          vel_advected[ijk + sp5_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm4_ed]);

                        const double vel_adv = interp_10(vel_advecting[ijk + sm4_ing],
                                                         vel_advecting[ijk + sm3_ing],
                                                         vel_advecting[ijk + sm2_ing],
                                                         vel_advecting[ijk + sm1_ing],
                                                         vel_advecting[ijk],
                                                         vel_advecting[ijk + sp1_ing],
                                                         vel_advecting[ijk + sp2_ing],
                                                         vel_advecting[ijk + sp3_ing],
                                                         vel_advecting[ijk + sp4_ing],
                                                         vel_advecting[ijk + sp5_ing]);

                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                    }
                }
            }
        }
        else{
            for(ssize_t i=imin;i<imax;i++){
                const ssize_t ishift = i*istride;
                for(ssize_t j=jmin;j<jmax;j++){
                    const ssize_t jshift = j*jstride;
                    for(ssize_t k=kmin;k<kmax;k++){
                        const ssize_t ijk = ishift + jshift + k;

                        //Upwind for positive velocity
                        const double phip = interp_weno11(vel_advected[ijk + sm5_ed],
                                                          vel_advected[ijk + sm4_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp5_ed]);

                        // Upwind for negative velocity
                        const double phim = interp_weno11(vel_advected[ijk + sp6_ed],
                                                          vel_advected[ijk + sp5_ed],
                                                          vel_advected[ijk + sp4_ed],
                                                          vel_advected[ijk + sp3_ed],
                                                          vel_advected[ijk + sp2_ed],
                                                          vel_advected[ijk + sp1_ed],
                                                          vel_advected[ijk],
                                                          vel_advected[ijk + sm1_ed],
                                                          vel_advected[ijk + sm2_ed],
                                                          vel_advected[ijk + sm3_ed],
                                                          vel_advected[ijk + sm4_ed]);

                        const double vel_adv = interp_10(vel_advecting[ijk + sm4_ing],
                                                         vel_advecting[ijk + sm3_ing],
                                                         vel_advecting[ijk + sm2_ing],
                                                         vel_advecting[ijk + sm1_ing],
                                                         vel_advecting[ijk],
                                                         vel_advecting[ijk + sp1_ing],
                                                         vel_advecting[ijk + sp2_ing],
                                                         vel_advecting[ijk + sp3_ing],
                                                         vel_advecting[ijk + sp4_ing],
                                                         vel_advecting[ijk + sp5_ing]);
                        flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                    }
                }
            }
        }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

    void fourth_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

        // Dynamically allocate flux array
        double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 1;
            const ssize_t jmin = 1;
            const ssize_t kmin = 1;

            const ssize_t imax = dims->nlg[0]-2;
            const ssize_t jmax = dims->nlg[1]-2;
            const ssize_t kmax = dims->nlg[2]-2;

            const ssize_t stencil[3] = {istride,jstride,1};
            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2 * sp1_ed ;
            const ssize_t sm1_ed = -sp1_ed ;

            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2 * sp1_ing;
            const ssize_t sm1_ing = -sp1_ing;

            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_4_pt(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                     interp_4_pt(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k];
                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_4_pt(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                     interp_4_pt(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed])) * rho0_half[k+1];
                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_4_pt(vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing]) *
                                     interp_4_pt(vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed]))* rho0[k];
                        }
                    }
                }
            }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

    void sixth_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

            // Dynamically allocate flux array
            double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 2;
            const ssize_t jmin = 2;
            const ssize_t kmin = 2;

            const ssize_t imax = dims->nlg[0]-3;
            const ssize_t jmax = dims->nlg[1]-3;
            const ssize_t kmax = dims->nlg[2]-3;

            const ssize_t stencil[3] = {istride,jstride,1};
            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2 * sp1_ed ;
            const ssize_t sp3_ed = 3 * sp1_ed ;
            const ssize_t sm1_ed = -sp1_ed ;
            const ssize_t sm2_ed = -2*sp1_ed;

            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2 * sp1_ing;
            const ssize_t sp3_ing = 3 * sp1_ing;
            const ssize_t sm1_ing = -sp1_ing;
            const ssize_t sm2_ing = -2*sp1_ing;

            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_6_pt(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                     interp_6_pt(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed])) * rho0_half[k];
                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_6_pt(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                     interp_6_pt(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed]))* rho0_half[k+1];
                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_6_pt(vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing]) *
                                     interp_6_pt(vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed]))* rho0[k];
                        }
                    }
                }
            }
        free(flux);
        return;
    }

    void eighth_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

            // Dynamically allocate flux array
            double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 3;
            const ssize_t jmin = 3;
            const ssize_t kmin = 3;

            const ssize_t imax = dims->nlg[0]-4;
            const ssize_t jmax = dims->nlg[1]-4;
            const ssize_t kmax = dims->nlg[2]-4;

            const ssize_t stencil[3] = {istride,jstride,1};
            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2 * sp1_ed ;
            const ssize_t sp3_ed = 3 * sp1_ed ;
            const ssize_t sp4_ed = 4 * sp1_ed ;
            const ssize_t sm1_ed = -sp1_ed ;
            const ssize_t sm2_ed = -2*sp1_ed;
            const ssize_t sm3_ed = -3*sp1_ed;

            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2 * sp1_ing;
            const ssize_t sp3_ing = 3 * sp1_ing;
            const ssize_t sp4_ing = 4 * sp1_ing;
            const ssize_t sm1_ing = -sp1_ing;
            const ssize_t sm2_ing = -2*sp1_ing;
            const ssize_t sm3_ing = -3*sp1_ing;

            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_8_pt(vel_advecting[ijk+sm3_ing],vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing],vel_advecting[ijk+sp4_ing]) *
                                     interp_8_pt(vel_advected[ijk+sm3_ed],vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed],vel_advected[ijk+sp4_ed])) * rho0_half[k];
                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_8_pt(vel_advecting[ijk+sm3_ing],vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing],vel_advecting[ijk+sp4_ing]) *
                                     interp_8_pt(vel_advected[ijk+sm3_ed],vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed],vel_advected[ijk+sp4_ed])) * rho0_half[k+1];
                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                        flux[ijk] = (interp_8_pt(vel_advecting[ijk+sm3_ing],vel_advecting[ijk+sm2_ing],vel_advecting[ijk+sm1_ing],vel_advecting[ijk],vel_advecting[ijk+sp1_ing],vel_advecting[ijk+sp2_ing],vel_advecting[ijk+sp3_ing],vel_advecting[ijk+sp4_ing]) *
                                     interp_8_pt(vel_advected[ijk+sm3_ed],vel_advected[ijk+sm2_ed],vel_advected[ijk+sm1_ed],vel_advected[ijk],vel_advected[ijk+sp1_ed],vel_advected[ijk+sp2_ed],vel_advected[ijk+sp3_ed],vel_advected[ijk+sp4_ed])) * rho0[k];
                        }
                    }
                }
            }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

    void weno_fifth_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

            // Dynamically allocate flux array
            double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 2;
            const ssize_t jmin = 2;
            const ssize_t kmin = 2;

            const ssize_t imax = dims->nlg[0]-3;
            const ssize_t jmax = dims->nlg[1]-3;
            const ssize_t kmax = dims->nlg[2]-3;

            const ssize_t stencil[3] = {istride,jstride,1};

            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2*sp1_ed ;
            const ssize_t sp3_ed = 3*sp1_ed ;
            const ssize_t sm1_ed = -sp1_ed ;
            const ssize_t sm2_ed = -2*sp1_ed ;

            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2*sp1_ing ;
            const ssize_t sm1_ing = -sp1_ing ;


            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed]);

                            const double vel_adv = interp_4_pt(vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed]);

                            const double vel_adv = interp_4_pt(vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno5(vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno5(vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed]);

                            const double vel_adv = interp_4_pt(vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                        }
                    }
                }
            }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }

    void weno_seventh_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

            // Dynamically allocate flux array
            double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 3;
            const ssize_t jmin = 3;
            const ssize_t kmin = 3;

            const ssize_t imax = dims->nlg[0]-4;
            const ssize_t jmax = dims->nlg[1]-4;
            const ssize_t kmax = dims->nlg[2]-4;

            const ssize_t stencil[3] = {istride,jstride,1};

            const ssize_t sm3_ed = -3*stencil[d_advecting];
            const ssize_t sm2_ed = -2*stencil[d_advecting];
            const ssize_t sm1_ed = -stencil[d_advecting];
            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2*stencil[d_advecting];
            const ssize_t sp3_ed = 3*stencil[d_advecting];
            const ssize_t sp4_ed = 4*stencil[d_advecting];

            const ssize_t sm2_ing = -2*stencil[d_advected];
            const ssize_t sm1_ing = -stencil[d_advected];
            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2*stencil[d_advected];
            const ssize_t sp3_ing = 3*stencil[d_advected];

            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp3_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed]);

                            const double vel_adv = interp_6_pt(vel_advecting[ijk + sm2_ing],
                                                            vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;
                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp3_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed]);

                            const double vel_adv = interp_6_pt(vel_advecting[ijk + sm2_ing],
                                                            vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];
                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;
                            //Upwind for positive velocity
                            const double phip = interp_weno7(vel_advected[ijk + sm3_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp3_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno7(vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed]);

                            const double vel_adv = interp_6_pt(vel_advecting[ijk + sm2_ing],
                                                            vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];
                        }
                    }
                }
            }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }


    void weno_ninth_order_m_pt(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
            double* restrict alpha0, double* restrict alpha0_half,
            double* restrict vel_advected, double* restrict vel_advecting,
            double* restrict tendency, ssize_t d_advected, ssize_t d_advecting){

            // Dynamically allocate flux array
            double *flux = (double *)malloc(sizeof(double)*dims->nlg[0] * dims->nlg[1] * dims->nlg[2]);

            const ssize_t istride = dims->nlg[1] * dims->nlg[2];
            const ssize_t jstride = dims->nlg[2];

            const ssize_t imin = 4;
            const ssize_t jmin = 4;
            const ssize_t kmin = 4;

            const ssize_t imax = dims->nlg[0]-5;
            const ssize_t jmax = dims->nlg[1]-5;
            const ssize_t kmax = dims->nlg[2]-5;

            const ssize_t stencil[3] = {istride,jstride,1};

            const ssize_t sm4_ed = -4*stencil[d_advecting];
            const ssize_t sm3_ed = -3*stencil[d_advecting];
            const ssize_t sm2_ed = -2*stencil[d_advecting];
            const ssize_t sm1_ed = -stencil[d_advecting];
            const ssize_t sp1_ed = stencil[d_advecting];
            const ssize_t sp2_ed = 2*stencil[d_advecting];
            const ssize_t sp3_ed = 3*stencil[d_advecting];
            const ssize_t sp4_ed = 4*stencil[d_advecting];
            const ssize_t sp5_ed = 5*stencil[d_advecting];

            const ssize_t sm3_ing = -3*stencil[d_advected];
            const ssize_t sm2_ing = -2*stencil[d_advected];
            const ssize_t sm1_ing = -stencil[d_advected];
            const ssize_t sp1_ing = stencil[d_advected];
            const ssize_t sp2_ing = 2*stencil[d_advected];
            const ssize_t sp3_ing = 3*stencil[d_advected];
            const ssize_t sp4_ing = 4*stencil[d_advected];

            if (d_advected != 2 && d_advecting !=2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                            vel_advected[ijk + sm3_ed],
                                                            vel_advected[ijk + sm2_ed],
                                                            vel_advected[ijk + sm1_ed],
                                                            vel_advected[ijk],
                                                            vel_advected[ijk + sp1_ed],
                                                            vel_advected[ijk + sp2_ed],
                                                            vel_advected[ijk + sp3_ed],
                                                            vel_advected[ijk + sp4_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                             vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm3_ed]);

                            const double vel_adv = interp_8_pt(vel_advecting[ijk + sm3_ing],
                                                            vel_advecting[ijk + sm2_ing],
                                                            vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing],
                                                            vel_advecting[ijk + sp4_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k] ;


                        }
                    }
                }
            }
            else if(d_advected == 2 && d_advecting == 2){
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                            vel_advected[ijk + sm3_ed],
                                                            vel_advected[ijk + sm2_ed],
                                                            vel_advected[ijk + sm1_ed],
                                                            vel_advected[ijk],
                                                            vel_advected[ijk + sp1_ed],
                                                            vel_advected[ijk + sp2_ed],
                                                            vel_advected[ijk + sp3_ed],
                                                            vel_advected[ijk + sp4_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                             vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm3_ed]);

                            const double vel_adv = interp_8_pt(vel_advecting[ijk+sm3_ing],
                                                            vel_advecting[ijk+sm2_ing],
                                                            vel_advecting[ijk+sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing],
                                                            vel_advecting[ijk + sp4_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0_half[k+1];

                        }
                    }
                }
            }
            else{
                for(ssize_t i=imin;i<imax;i++){
                    const ssize_t ishift = i*istride;
                    for(ssize_t j=jmin;j<jmax;j++){
                        const ssize_t jshift = j*jstride;
                        for(ssize_t k=kmin;k<kmax;k++){
                            const ssize_t ijk = ishift + jshift + k;

                            //Upwind for positive velocity
                            const double phip = interp_weno9(vel_advected[ijk + sm4_ed],
                                                            vel_advected[ijk + sm3_ed],
                                                            vel_advected[ijk + sm2_ed],
                                                            vel_advected[ijk + sm1_ed],
                                                            vel_advected[ijk],
                                                            vel_advected[ijk + sp1_ed],
                                                            vel_advected[ijk + sp2_ed],
                                                            vel_advected[ijk + sp3_ed],
                                                            vel_advected[ijk + sp4_ed]);

                            // Upwind for negative velocity
                            const double phim = interp_weno9(vel_advected[ijk + sp5_ed],
                                                             vel_advected[ijk + sp4_ed],
                                                             vel_advected[ijk + sp3_ed],
                                                             vel_advected[ijk + sp2_ed],
                                                             vel_advected[ijk + sp1_ed],
                                                             vel_advected[ijk],
                                                             vel_advected[ijk + sm1_ed],
                                                             vel_advected[ijk + sm2_ed],
                                                             vel_advected[ijk + sm3_ed]);

                            const double vel_adv = interp_8_pt(vel_advecting[ijk + sm3_ing],
                                                            vel_advecting[ijk + sm2_ing],
                                                            vel_advecting[ijk + sm1_ing],
                                                            vel_advecting[ijk],
                                                            vel_advecting[ijk + sp1_ing],
                                                            vel_advecting[ijk + sp2_ing],
                                                            vel_advecting[ijk + sp3_ing],
                                                            vel_advecting[ijk + sp4_ing]);

                            flux[ijk] = 0.5 * ((vel_adv+fabs(vel_adv))*phip + (vel_adv-fabs(vel_adv))*phim)*rho0[k];

                        }
                    }
                }
            }
        momentum_flux_divergence(dims, alpha0, alpha0_half, flux,
                                tendency, d_advected, d_advecting);
        free(flux);
        return;
    }


void compute_advective_tendencies_m(struct DimStruct *dims, double* restrict rho0, double* restrict rho0_half,
                                double* restrict alpha0, double* restrict alpha0_half,
                                double* restrict vel_advected, double* restrict vel_advecting,
                                double* restrict tendency, ssize_t d_advected, ssize_t d_advecting, int scheme){



    switch(scheme){
        case 2:
            second_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 3:
            weno_third_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 4:
            fourth_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 5:
            weno_fifth_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 6:
           sixth_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 7:
            weno_seventh_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 8:
            eighth_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 9:
            weno_ninth_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 11:
            weno_eleventh_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 14:
            // This is an application of fourth order Wicker-Skamarock to momentum but using a lower order interpolation
            // for advecting velocity.
            fourth_order_ws_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 16:
            // This is an application of sixth order Wicker-Skamarock to momentum but using a lower order interpolation
            // for advecting velocity.
            sixth_order_ws_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 24:
            fourth_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 25:
            weno_fifth_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 26:
            sixth_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 27:
            weno_seventh_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 28:
            eighth_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 29:
            weno_ninth_order_m_pt(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 32:
            second_order_ml_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 34:
            fourth_order_ml_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        case 36:
            sixth_order_ml_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
        default:
            // Default to second order scheme.
            second_order_m(dims, rho0, rho0_half, alpha0, alpha0_half, vel_advected, vel_advecting,
                tendency, d_advected, d_advecting);
            break;
    };
}
