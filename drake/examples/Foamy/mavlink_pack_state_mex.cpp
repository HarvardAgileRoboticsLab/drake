#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>
#include <math.h>

#define DEBUG_PRINTING

#define PI 3.141592653589793
#define RE 6371000.0 //Earth radius in meters
#define BASELAT 42.467332 //Base station latitude in degrees
#define COSLAT 0.737662414272813 //cosine of base station latitude
#define BASELON -71.414699 //Base station longitude in degrees
#define BASEALT 413.0 //Base station longitude in degrees

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
    if(len_y < 15) {
        mexErrMsgTxt("Second input must be a double array with 15 elements.");
        return;
    }
    
    double xenu = x[0];
    double yenu = x[1];
    double zenu = x[2];
    
    double lat = round(10000000.0*(BASELAT + (180.0/PI)*(yenu/RE))); //latitude in degrees*10^-7 (linearized)
    double lon = round(10000000.0*(BASELON + (180.0/PI)*(xenu/(RE*COSLAT)))); //longitude in degrees*10^-7 (linearized)
    double alt = round(1000.0*(BASEALT + zenu)); //altitude in millimeters
    
    float attitude_quaternion[4];
    attitude_quaternion[0] = (float)x[3];
    attitude_quaternion[1] = (float)x[4];
    attitude_quaternion[2] = (float)x[5];
    attitude_quaternion[3] = (float)x[6];
    
    double xdotenu = x[7];
    double ydotenu = x[8];
    double zdotenu = x[9];
    
    double vlat = round(100.0*ydotenu); //linear velocity in cm/sec
    double vlon = round(100.0*xdotenu); //linear velocity in cm/sec
    double valt = round(100.0*zdotenu); //linear velocity in cm/sec
    
    double true_airspeed = round(sqrt(vlat*vlat + vlon+vlon + valt*valt));
    double ind_airspeed = true_airspeed;
    
    float wx = (float)x[10];
    float wy = (float)x[11];
    float wz = (float)x[12];
    
    double xacc = round((1000/9.81)*y[6]); //convert to g*10^-3
    double yacc = round((1000/9.81)*y[7]); //convert to g*10^-3
    double zacc = round((1000/9.81)*y[8]); //convert to g*10^-3
    
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