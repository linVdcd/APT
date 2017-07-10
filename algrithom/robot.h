/*-----------------------------------
	林明安 2016.6.1
	两种机器人音效，1、基于环形调制；2、基于相位置零。
*/
#ifndef ROBOT_H
#define ROBOT_H
#include "Platform.h"

class RobotSound{

	
public:
	
	
	void Robot(const short *x, int Lx,short *y, const short type, const int fs,const double RMfre,const int outwin);
	void WhisperSound(const short *x,int Lx, short **y, int &Ly,const int p);
	RobotSound();
	~RobotSound();
};
#endif