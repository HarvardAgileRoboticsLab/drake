#include "mex.h"
#include <mavlink.h>

#define DEBUG_PRINTING

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    //Do some checks
    if(nrhs < 1) {
        mexErrMsgTxt("Not enough input arguments.");
        return;
    }
    if(nrhs > 1) {
        mexErrMsgTxt("Too many input arguments.");
        return;
    }
    if(nlhs > 2) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
    if(!(mxGetClassID(prhs[0]) == mxUINT8_CLASS)) {
        mexErrMsgTxt("Input must be a byte (uint8) array.");
        return;
    }
    
    uint8_t *p = (uint8_t *)mxGetData(prhs[0]);
    mwSize len = mxGetNumberOfElements(prhs[0]);
    
    mavlink_message_t rxbuf;
    mavlink_status_t sbuf;

    mavlink_message_t msg;
    mavlink_status_t sout;

    uint8_t received = 0;
    for(int k = 0; k < len; k++) {
        received = mavlink_frame_char_buffer(&rxbuf, &sbuf, p[k], &msg, &sout);
    }
    
    #ifdef DEBUG_PRINTING
    if(received == 0) {
        mexPrintf("Not received.\n");
    }
    else if(received == 2) {
        mexPrintf("CRC error.\n");
    }
    mexPrintf("Message ID: %d\n", msg.msgid);
    #endif
    
    if(nlhs > 0 && received == 1) {
        
        if(msg.msgid == MAVLINK_MSG_ID_HIL_ACTUATOR_CONTROLS) {
            plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
            double *u = mxGetPr(plhs[0]);
            float controls[16];
            mavlink_msg_hil_actuator_controls_get_controls(&msg, controls);
            //These mappings come from plane_sitl.main.mix
            u[0] = (double)controls[4]; //throttle
            u[1] = (double)controls[6]; //aileron (roll)
            u[2] = (double)controls[7]; //elevator (pitch)
            u[3] = (double)controls[2]; //rudder (yaw)
        }
        else {
            for(int k = 0; k < nlhs; k++) {
                plhs[k] = mxCreateDoubleMatrix(0, 0, mxREAL);
            }
            return;
        }
        
        if(nlhs > 1) {
            plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
            double *t = mxGetPr(plhs[1]);
            t[0] = (double)mavlink_msg_hil_controls_get_time_usec(&msg);
        }

    }
    else {
        for(int k = 0; k < nlhs; k++) {
            plhs[k] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
    }
    
    return;
}