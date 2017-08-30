#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>
#include <math.h>

#define SQRT2 1.414213562373095

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
    if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])) {
        mexErrMsgTxt("Inputs must be double arrays.");
        return;
    }
    
    double *x = mxGetPr(prhs[0]);
    mwSize len_x = mxGetNumberOfElements(prhs[0]);
    if(len_x < 13) {
        mexErrMsgTxt("First input must be a double array with 13 elements.");
        return;
    }
    
    double *y = mxGetPr(prhs[1]);
    mwSize len_y = mxGetNumberOfElements(prhs[1]);
    if(len_y < 16) {
        mexErrMsgTxt("Second input must be a double array with 17 elements.");
        return;
    }
    
    double lat = round(10000000.0*y[0]); //latitude in degrees*10^-7
    double lon = round(10000000.0*y[1]); //longitude in degrees*10^-7
    double alt = round(1000.0*y[2]); //altitude in millimeters
    
    //We're using ENU in our simulator but PX4 uses NED
    //So we have to rotate our quaternion
    float attitude_quaternion[4];
    attitude_quaternion[0] = (float)((1/SQRT2)*(x[4]+x[5]));
    attitude_quaternion[1] = (float)((-1/SQRT2)*(x[3]+x[6]));
    attitude_quaternion[2] = (float)((-1/SQRT2)*(x[3]-x[6]));
    attitude_quaternion[3] = (float)((-1/SQRT2)*(x[5]-x[4]));
    
    double vlat = round(100.0*y[3]); //linear velocity in cm/sec
    double vlon = round(100.0*y[4]); //linear velocity in cm/sec
    double valt = round(-100.0*y[5]); //linear velocity in cm/sec
    
    double true_airspeed = round(sqrt(vlat*vlat + vlon+vlon + valt*valt));
    double ind_airspeed = true_airspeed;
    
    float wx = (float)x[10];
    float wy = (float)x[11];
    float wz = (float)x[12];
    
    double xacc = 0; //round((1000/9.81)*y[6]); //convert to g*10^-3
    double yacc = 0; //round((1000/9.81)*y[7]); //convert to g*10^-3
    double zacc = 0; //round((1000/9.81)*y[8]); //convert to g*10^-3
    
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
    
    timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t time_usec = 1000000*tv.tv_sec + tv.tv_usec;
    
    len = mavlink_msg_hil_state_quaternion_pack(0x01, 0xc8, &msg, time_usec,
                                                attitude_quaternion, //float*
                                                wx, wy, wz, //float (rad/s)
                                                (int32_t)lat, (int32_t)lon, //int32 (deg*10^-7)
                                                (int32_t)alt, //int32 (millimeters)
                                                (int16_t)vlat, (int16_t)vlon, (int16_t)valt, //int16 (cm/sec)
                                                (int16_t)ind_airspeed, //int16 (cm/sec)
                                                (int16_t)true_airspeed, //int16 (cm/sec)
                                                (int16_t)xacc, (int16_t)yacc, (int16_t)zacc); //int16 (g*10^-3)
    
    len = mavlink_msg_to_send_buffer(buf, &msg);
    
    plhs[0] = mxCreateNumericMatrix(len, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *s = (uint8_t *)mxGetData(plhs[0]);
    
    for(int k = 0; k < len; k++) {
        s[k] = buf[k];
    }
    
    return;
}