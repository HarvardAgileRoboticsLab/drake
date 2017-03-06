#include "cart_pole.h"
#include <math.h>

#define pi M_PI

const double CartPole::mc=10;
const double CartPole::mp=1;
const double CartPole::l=0.5;
const double CartPole::g=9.81;

CartPole::CartPole(double& mydt, unsigned int& myN)
{
    stateNb=4;
    commandNb=1;
    dt = mydt;
    N = myN;
    Id.setIdentity();

    fxList.resize(N);
    fuList.resize(N);
    
    fxx_new[0].setZero();
    fxx_new[1].setZero();
    fxx_new[2].setZero();
    fxx_new[3].setZero();
    fuu_new[0].setZero();
    fux_new[0].setZero();
    fxu_new[0].setZero();

    lowerCommandBounds << -50.0;
    upperCommandBounds << 50.0;

    H.setZero();
    C.setZero();
    G.setZero();
    Bu.setZero();
    A1.setZero();
    A2.setZero();
    A3.setZero();
    A4.setZero();
    B1.setZero();
    B2.setZero();
    B3.setZero();
    B4.setZero();
    IdentityMat.setIdentity();

    Xp1.setZero();
    Xp2.setZero();
    Xp3.setZero();
    Xp4.setZero();

    Xm1.setZero();
    Xm2.setZero();
    Xm3.setZero();
    Xm4.setZero();

    debugging_print = 0;
}

stateVec_t CartPole::cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U)
{
    H << mc + mp, mp*l*cos(X(1,0)),
         mp*l*cos(X(1,0)), mp*pow(l,2.0);
    C << 0, -mp*X(3,0)*l*sin(X(1,0)),
         0, 0;
    G << 0,
         mp*g*l*sin(X(1,0));
    Bu << 1,
         0;     
    stateVec_half_t velocity;
    velocity << X(2),
                X(3);
    stateVec_half_t accel = - H.inverse() * (C*velocity + G - Bu*U);

    stateVec_t X_new;
    X_new << velocity(0),
             velocity(1),
             accel(0),
             accel(1);

    return X_new;
}

void CartPole::cart_pole_dyn_cst(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, 
                                stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    // costFunctionCartPole
    if(debugging_print) TRACE_CART_POLE("initialize dimensions\n");
    int N = xList.size();//[TODO: to be checked]
    int n = xList[0].rows();
    int m = uList[0].rows();
    CostFunctionCartPole *costFunctionCartPole = &costFunction_cart_pole;

    c = 0;
    stateMat_t AA;//not used
    stateVec_t BB;//not used
    AA.setZero();
    BB.setZero();

    if(debugging_print) TRACE_CART_POLE("compute cost function\n");

    commandMat_t c_mat_to_scalar;

    stateMatTab_t A_temp;//not used
    stateR_commandC_tab_t B_temp;//not used
    A_temp.resize(N);
    B_temp.resize(N);
    
    if(nargout == 2){
        const int nargout_update1 = 3;        
        for(unsigned int k=0;k<N;k++){
            if(isNan(uList[k])){
                if(debugging_print) TRACE_CART_POLE("before the update1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose()) * costFunctionCartPole->getQf() * (xList[k] - xgoal);
                c += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_CART_POLE("after the update1\n");
            }else{
                if(debugging_print) TRACE_CART_POLE("before the update2\n");
                FList[k] = update(nargout_update1, dt, xList[k], uList[k], AA, BB);
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_CART_POLE("after the update2\n");
                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunctionCartPole->getR()*uList[k];
                c += c_mat_to_scalar(0,0);
            }
        }
    }else{
        const int nargout_update2 = 3;
        for(unsigned int k=0;k<N;k++){
            if(isNan(uList[k])){
                if(debugging_print) TRACE_CART_POLE("before the update3\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQf()*(xList[k] - xgoal);
                c += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_CART_POLE("after the update3\n");
            }else{
                if(debugging_print) TRACE_CART_POLE("before the update4\n");
                FList[k] = update(nargout_update2, dt, xList[k], uList[k], AA, BB);//assume three outputs, code needs to be optimized
                if(debugging_print) TRACE_CART_POLE("before the update4-1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_CART_POLE("after the update4\n");

                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunctionCartPole->getR()*uList[k];
                c += c_mat_to_scalar(0,0); // TODO: to be checked
                if(debugging_print) TRACE_CART_POLE("after the update5\n");
                
                A_temp[k] = AA;
                B_temp[k] = BB;
            }
        }

        stateVec_t cx_temp;
        
        if(debugging_print) TRACE_CART_POLE("compute dynamics and cost derivative\n");

        for(unsigned int k=0;k<N-1;k++){
            fxList[k] = A_temp[k];
            fuList[k] = B_temp[k];
            
            cx_temp << xList[k](0,0)-xgoal(0), xList[k](1,0)-xgoal(1), xList[k](2,0)-xgoal(2), xList[k](3,0)-xgoal(3);
            cx[k] = costFunctionCartPole->getQ()*cx_temp;
            cu[k] = costFunctionCartPole->getR()*uList[k];
            cxx[k] = costFunctionCartPole->getQ();
            cux[k].setZero();
            cuu[k] = costFunctionCartPole->getR();
            
            // if(k == 49) {
            //     std::cout << "fxList[49]: " << fxList[k] << std::endl; 
            //     std::cout << "fuList[49]: " << fuList[k] << std::endl;   
            //     std::cout << "cx[49]: " << cx[k] << std::endl;
            //     std::cout << "cu[49]: " << cu[k] << std::endl;
            //     std::cout << "cxx[49]: " << cxx[k] << std::endl;
            // }

        }
        if(debugging_print) TRACE_CART_POLE("update the final value of cost derivative \n");

        cx[N-1] = costFunctionCartPole->getQf()*(xList[N-1]-xgoal);//[TODO: double check whether there is - xgoal]
        cu[N-1] = costFunctionCartPole->getR()*uList[N-1];
        cxx[N-1] = costFunctionCartPole->getQf();
        cux[N-1].setZero();
        cuu[N-1] = costFunctionCartPole->getR();

        if(debugging_print) TRACE_CART_POLE("set unused matrices to zero \n");

        // the following matrices and scalars are set to Zero instead of empty, not supported by Eigen.
        //fxx, fxu, fuu are not defined since never used
        for(unsigned int k=0;k<N;k++){
            FList[k].setZero();
        }    
        c = 0;
    }
    if(debugging_print) TRACE_CART_POLE("finish cart_pole_dyn_cst\n");
}

void CartPole::cart_pole_dyn_cst_short(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, const stateVec_t& xgoal, stateVec_t& xList_next, double& c){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    // costFunctionCartPole
    if(debugging_print) TRACE_CART_POLE("initialize dimensions\n");
    int N = xList_curr.cols();//[TODO: to be checked, should be 1 in this case, the for loop only runs once, the code needs to be optimized]
    int n = xList_curr.rows();
    int m = uList_curr.rows();
    CostFunctionCartPole *costFunctionCartPole_short = &costFunction_cart_pole;

    c = 0;
    stateMat_t AA;
    stateVec_t BB;
    AA.setZero();
    BB.setZero();

    if(debugging_print) TRACE_CART_POLE("compute cost function\n");

    commandMat_t c_mat_to_scalar;

    stateMatTab_t A_temp;
    stateR_commandC_tab_t B_temp;
    A_temp.resize(N);
    B_temp.resize(N);
    xList_next.setZero();

    const int nargout_update1 = 1;
    for(unsigned int k=0;k<N;k++){
        if(isNan(uList_curr)){ 
            if(debugging_print) TRACE_CART_POLE("before the update1\n");
            c_mat_to_scalar = 0.5*(xList_curr.transpose() - xgoal.transpose()) * costFunctionCartPole_short->getQf() * (xList_curr - xgoal);
            // std::cout << "size of N: " << N << std::endl;
            //std::cout << "x_final_state: " << xList_curr.transpose() << std::endl;
            c += c_mat_to_scalar(0,0);
            //std::cout << "c: " << c << std::endl;
            if(debugging_print) TRACE_CART_POLE("after the update1\n");
        }else{
            if(debugging_print) TRACE_CART_POLE("before the update2\n");
            xList_next = update(nargout_update1, dt, xList_curr, uList_curr, AA, BB);
            c_mat_to_scalar = 0.5*(xList_curr.transpose() - xgoal.transpose())*costFunctionCartPole_short->getQ()*(xList_curr - xgoal);
            //if(k<5) std::cout << "c_mat_to_scalar1: " << c_mat_to_scalar << std::endl;
            if(debugging_print) TRACE_CART_POLE("after the update2\n");
            c_mat_to_scalar += 0.5*uList_curr.transpose()*costFunctionCartPole_short->getR()*uList_curr;
            //if(k<5) std::cout << "c_mat_to_scalar2: " << c_mat_to_scalar << std::endl;
            c += c_mat_to_scalar(0,0);
        }
    }
    if(debugging_print) TRACE_CART_POLE("finish cart_pole_dyn_cst\n");
}

stateVec_t CartPole::update(const int& nargout, const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B){
    // 4th-order Runge-Kutta step
    if(debugging_print) TRACE_CART_POLE("update: 4th-order Runge-Kutta step\n");
    Xdot1 = cart_pole_dynamics(X, U);
    Xdot2 = cart_pole_dynamics(X + 0.5*dt*Xdot1, U);
    Xdot3 = cart_pole_dynamics(X + 0.5*dt*Xdot2, U);
    Xdot4 = cart_pole_dynamics(X + dt*Xdot3, U);
    stateVec_t X_new;
    X_new = X + (dt/6)*(Xdot1 + 2*Xdot2 + 2*Xdot3 + Xdot4);
    
    if(debugging_print) TRACE_CART_POLE("update: X_new\n");

    if(nargout > 1){
        int n = X.size();
        int m = U.size();

        double delta = 1e-7;
        stateMat_t Dx;
        commandMat_t Du;
        Dx.setIdentity();
        Dx = delta*Dx;
        Du.setIdentity();
        Du = delta*Du;

        for(unsigned int i=0;i<n;i++){
            Xp1 = cart_pole_dynamics(X+Dx.col(i),U);
            Xm1 = cart_pole_dynamics(X-Dx.col(i),U);
            A1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = cart_pole_dynamics(X+0.5*dt*Xdot1+Dx.col(i),U);
            Xm2 = cart_pole_dynamics(X+0.5*dt*Xdot1-Dx.col(i),U);
            A2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = cart_pole_dynamics(X+0.5*dt*Xdot2+Dx.col(i),U);
            Xm3 = cart_pole_dynamics(X+0.5*dt*Xdot2-Dx.col(i),U);
            A3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = cart_pole_dynamics(X+0.5*dt*Xdot3+Dx.col(i),U);
            Xm4 = cart_pole_dynamics(X+0.5*dt*Xdot3-Dx.col(i),U);
            A4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        for(unsigned int i=0;i<m;i++){
            Xp1 = cart_pole_dynamics(X,U+Du.col(i));
            Xm1 = cart_pole_dynamics(X,U-Du.col(i));
            B1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = cart_pole_dynamics(X+0.5*dt*Xdot1,U+Du.col(i));
            Xm2 = cart_pole_dynamics(X+0.5*dt*Xdot1,U-Du.col(i));
            B2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = cart_pole_dynamics(X+0.5*dt*Xdot2,U+Du.col(i));
            Xm3 = cart_pole_dynamics(X+0.5*dt*Xdot2,U-Du.col(i));
            B3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = cart_pole_dynamics(X+0.5*dt*Xdot3,U+Du.col(i));
            Xm4 = cart_pole_dynamics(X+0.5*dt*Xdot3,U-Du.col(i));
            B4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        A = (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)*(IdentityMat + A2 * dt/3)*(IdentityMat + A1 * dt/6);
        B = B4 * dt/6 + (IdentityMat + A4 * dt/6) * B3 * dt/3 + (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)* B2 * dt/3 + (IdentityMat + (dt/6)*A4)*(IdentityMat + (dt/3)*A3)*(IdentityMat + (dt/3)*A2)*(dt/6)*B1;
    }
    if(debugging_print) TRACE_CART_POLE("update: X_new\n");
    return X_new;
}

void CartPole::grad(const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B){
    int n = X.size();
    int m = U.size();

    double delta = 1e-7;
    stateMat_t Dx;
    commandMat_t Du;
    Dx.setIdentity();
    Dx = delta*Dx;
    Du.setIdentity();
    Du = delta*Du;

    stateMat_t AA;
    stateVec_t BB;
    AA.setZero();
    BB.setZero();

    int nargout = 1;
    for(unsigned int i=0;i<n;i++){
        Xp = update(nargout, dt, X+Dx.col(i), U, AA, BB);
        Xm = update(nargout, dt, X-Dx.col(i), U, AA, BB);
        A.col(i) = (Xp - Xm)/(2*delta);
    }

    for(unsigned int i=0;i<m;i++){
        Xp = update(nargout, dt, X, U+Du.col(i), AA, BB);
        Xm = update(nargout, dt, X, U-Du.col(i), AA, BB);
        B.col(i) = (Xp - Xm)/(2*delta);
    }
}

void CartPole::hessian(const double& dt, const stateVec_t& X, const commandVec_t& U){
    int n = X.size();
    int m = U.size();

    double delta = 1e-5;
    stateMat_t Dx;
    commandMat_t Du;
    Dx.setIdentity();
    Dx = delta*Dx;
    Du.setIdentity();
    Du = delta*Du;

    stateMat_t Ap;
    Ap.setZero();
    stateMat_t Am;
    Am.setZero();
    stateVec_t B;
    B.setZero();

    for(unsigned int i=0;i<n;i++){
        fxx_new[i].setZero();
        fxu_new[i].setZero();
        fuu_new[i].setZero();
    }

    for(unsigned int i=0;i<n;i++){
        grad(dt, X+Dx.col(i), U, Ap, B);
        grad(dt, X-Dx.col(i), U, Am, B);
        fxx_new[i] = (Ap - Am)/(2*delta);
    }

    stateVec_t Bp;
    Bp.setZero();
    stateVec_t Bm;
    Bm.setZero();

    for(unsigned int j=0;j<m;j++){
        grad(dt, X, U+Du.col(j), Ap, Bp);
        grad(dt, X, U-Du.col(j), Am, Bm);
        fxu_new[j] = (Ap - Am)/(2*delta);
        fuu_new[j] = (Bp - Bm)/(2*delta);
    }
}