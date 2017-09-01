#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>
#include <math.h>

#define SQRT2 1.414213562373095
#define BASEALT 413.0

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
    if(nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
        return;
    }
    if(!mxIsDouble(prhs[0])) {
        mexErrMsgTxt("Input must be a double array.");
        return;
    }
    
    double *p = mxGetPr(prhs[0]);
    mwSize len_in = mxGetNumberOfElements(prhs[0]);
    
    if(len_in < 13) {
        mexErrMsgTxt("Input must be a double array with 13 elements.");
        return;
    }
    
    //Position in NED frame in meters
    double x = p[1];
    double y = p[0];
    double z = -p[2];
    
    //We're using ENU in our simulator but PX4 uses NED
    //So we have to rotate our quaternion
    double q0 = ((1/SQRT2)*(p[4]+p[5]));
    double q1 = ((-1/SQRT2)*(p[3]+p[6]));
    double q2 = ((-1/SQRT2)*(p[3]-p[6]));
    double q3 = ((-1/SQRT2)*(p[5]-p[4]));
    
    //Convert quaternion to Euler angles
    double roll = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1*q1 + q2*q2));
    double pitch = asin(2*(q0*q2 - q3*q1));
    double yaw = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2*q2 + q3*q3));
    
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
    
    timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t time_usec = 1000000*tv.tv_sec + tv.tv_usec;
    
    len = mavlink_msg_vision_position_estimate_pack(0x01, 0xc8, &msg, time_usec,
                                                    (float)x, (float)y, (float)z,
                                                    (float)roll, (float)pitch, (float)yaw);
    
    len = mavlink_msg_to_send_buffer(buf, &msg);
    
    plhs[0] = mxCreateNumericMatrix(len, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *s = (uint8_t *)mxGetData(plhs[0]);
    
    for(int k = 0; k < len; k++) {
        s[k] = buf[k];
    }
    
    return;
}