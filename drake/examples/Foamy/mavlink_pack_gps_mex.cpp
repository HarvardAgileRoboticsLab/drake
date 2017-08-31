#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>
#include <math.h>

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
    
    double *y = mxGetPr(prhs[0]);
    mwSize len_in = mxGetNumberOfElements(prhs[0]);
    
    if(len_in < 17) {
        mexErrMsgTxt("Input must be a double array with 17 elements.");
        return;
    }
    
    double lat = round(10000000.0*y[0]); //latitude in degrees*10^-7
    double lon = round(10000000.0*y[1]); //longitude in degrees*10^-7
    double alt = round(1000.0*y[2]); //altitude in millimeters
    
    double vn = round(100.0*y[3]); //velocity in cm/s in north direction
    double ve = round(100.0*y[4]); //velocity in cm/s in east direction
    double vd = round(100.0*y[5]); //velocity in cm/s in down direction
    
    double vel = round(100.0*sqrt(y[3]*y[3] + y[4]*y[4] + y[5]*y[5])); //velocity magnitude in cm/s
    
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
    
    timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t time_usec = 1000000*tv.tv_sec + tv.tv_usec;
    
    len = mavlink_msg_hil_gps_pack(0x01, 0xc8, &msg, time_usec, 3,
                                   (int32_t)lat, (int32_t)lon, (int32_t)alt,
                                   100, 100, (int16_t)vel,
                                   (int16_t)vn, (int16_t)ve, (int16_t)vd,
                                   65535, 10);
    
    len = mavlink_msg_to_send_buffer(buf, &msg);
    
    plhs[0] = mxCreateNumericMatrix(len, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *s = (uint8_t *)mxGetData(plhs[0]);
    
    for(int k = 0; k < len; k++) {
        s[k] = buf[k];
    }
    
    return;
}