function [b] = getFrictionTorqueVector(obj, qd)

bv1_positive=obj.bv1_positive; bv2_positive=obj.bv2_positive; bv3_positive=obj.bv3_positive; bv4_positive=obj.bv4_positive; bv5_positive=obj.bv5_positive; bv6_positive=obj.bv6_positive; bv7_positive=obj.bv7_positive;
bv1_negative=obj.bv1_negative; bv2_negative=obj.bv2_negative; bv3_negative=obj.bv3_negative; bv4_negative=obj.bv4_negative; bv5_negative=obj.bv5_negative; bv6_negative=obj.bv6_negative; bv7_negative=obj.bv7_negative;
bc1_positive=obj.bc1_positive; bc2_positive=obj.bc2_positive; bc3_positive=obj.bc3_positive; bc4_positive=obj.bc4_positive; bc5_positive=obj.bc5_positive; bc6_positive=obj.bc6_positive; bc7_positive=obj.bc7_positive;
bc1_negative=obj.bc1_negative; bc2_negative=obj.bc2_negative; bc3_negative=obj.bc3_negative; bc4_negative=obj.bc4_negative; bc5_negative=obj.bc5_negative; bc6_negative=obj.bc6_negative; bc7_negative=obj.bc7_negative;


% joint angular velocity data
dq1 = qd(1);
dq2 = qd(2);
dq3 = qd(3);
dq4 = qd(4);
dq5 = qd(5);
dq6 = qd(6);
dq7 = qd(7);

b1 = (bv1_positive + bv1_negative)/2 * dq(1);
b2 = (bv2_positive + bv2_negative)/2 * dq(2);
b3 = (bv3_positive + bv3_negative)/2 * dq(3);
b4 = (bv4_positive + bv4_negative)/2 * dq(4);
b5 = (bv5_positive + bv5_negative)/2 * dq(5);
b6 = (bv6_positive + bv6_negative)/2 * dq(6);
b7 = (bv7_positive + bv7_negative)/2 * dq(7);

b = [b1;b2;b3;b4;b5;b6;b7];
end