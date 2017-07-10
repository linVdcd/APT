#pragma once
#ifndef VIBRATO_H
#define VIBRATO_H

void VibratoEffect(const short *x, const int Lx, const int fs, short *y,const double modfreq, const double width);

#endif // !VIBRATO_H