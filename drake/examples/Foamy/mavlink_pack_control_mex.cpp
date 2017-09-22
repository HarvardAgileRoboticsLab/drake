#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //Do some checks
    if(nrhs < 2) {
        mexErrMsgTxt("Not enough input arguments.");
        return;
    }
    if(nrhs > 2) {
        mexErrMsgTxt("Too many input arguments.");
        return;
    }
    if(nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
    if(!mxIsDouble(prhs[0])) {
        mexErrMsgTxt("Input must be a double array.");
        return;
    }
    
    double *u = mxGetPr(prhs[0]);
    double *t = mxGetPr(prhs[1]);
    mwSize len_in = mxGetNumberOfElements(prhs[0]);
    
    if(len_in < 4) {
        mexErrMsgTxt("Input must be a double array with 4 elements.");
        return;
    }
    
    float thr = (float)u[0];
    float ail = (float)u[1];
    float ele = (float)u[2];
    float rud = (float)u[3];
    
    float controls[16] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    controls[4] = thr;
    controls[6] = ail;
    controls[5] = -ail;
    controls[7] = ele;
    controls[2] = rud;
    
    double time_usec = 1000000.0*t[0];

    uint8_t mode = 129;
    uint64_t flag = 0;
    
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
        
    len = mavlink_msg_hil_actuator_controls_pack(0x01, 0xc8, &msg, (uint64_t)time_usec,
                                                 controls, mode, flag);
    
    len = mavlink_msg_to_send_buffer(buf, &msg);
    
    plhs[0] = mxCreateNumericMatrix(len, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *s = (uint8_t *)mxGetData(plhs[0]);
    
    for(int k = 0; k < len; k++) {
        s[k] = buf[k];
    }
    
    return;
}