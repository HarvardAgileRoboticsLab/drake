#include "mex.h"
#include <mavlink.h>
#include <sys/time.h>

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
    
    if(len_in < 17) {
        mexErrMsgTxt("Input must be a double array with 17 elements.");
        return;
    }
    
    float xacc = (float)p[6];
    float yacc = (float)p[7];
    float zacc = (float)p[8];
    
    float xgyro = (float)p[9];
    float ygyro = (float)p[10];
    float zgyro = (float)p[11];
    
    float xmag = (float)p[12];
    float ymag = (float)p[13];
    float zmag = (float)p[14];
    
    float abs_pressure = (float)p[15];
    float diff_pressure = (float)p[16];
    float pressure_alt = (float)p[2];
    float temperature = 15.0;
    
    uint32_t fields_updated = 0x0fff;
    
    mavlink_message_t msg;
    uint8_t buf[MAVLINK_MAX_PACKET_LEN];
    uint16_t len;
    
    timeval tv;
    gettimeofday(&tv, NULL);
    uint64_t time_usec = 1000000*tv.tv_sec + tv.tv_usec;
    
    len = mavlink_msg_hil_sensor_pack(0x01, 0xc8, &msg, time_usec,
                                      xacc, yacc, zacc,
                                      xgyro, ygyro, zgyro,
                                      xmag, ymag, zmag,
                                      abs_pressure, diff_pressure,
                                      pressure_alt, temperature,
                                      fields_updated);
    
    len = mavlink_msg_to_send_buffer(buf, &msg);
    
    plhs[0] = mxCreateNumericMatrix(len, 1, mxUINT8_CLASS, mxREAL);
    uint8_t *s = (uint8_t *)mxGetData(plhs[0]);
    
    for(int k = 0; k < len; k++) {
        s[k] = buf[k];
    }
    
    return;
}