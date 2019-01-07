#include <math.h>
#include <iostream>

using namespace std;

typedef unsigned int uint;

const uint MAX_POWER = sizeof(uint)*8;

unsigned int bit_reverse(register uint x){
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return((x >> 16) | (x << 16));
}

// complex numbers
struct cfloat{
	float re, im;

	cfloat(float real = 0.0f, float imag = 0.0f): re(real), im(imag) {}

	cfloat operator+(const cfloat &c) const { return cfloat(re + c.re, im + c.im);	}
	cfloat operator-(const cfloat &c) const { return cfloat(re - c.re, im - c.im);	}
	cfloat operator*(const cfloat &c) const { return cfloat(re*c.re - im*c.im, re*c.im + im*c.re);	}
	friend cfloat operator*(float k, const cfloat &c){ return cfloat(k*c.re, k*c.im); }
};

// output for complex numbers
void show(const cfloat &c){ cout << "(" << c.re << " i " << c.im << ")" << endl; }


cfloat *twiddle_lut = 0;
uint twiddle_lut_power;

const float twiddle_epsilon = 1e-6f;

void initTwiddleLUT(uint power){
	if (twiddle_lut != 0)
		delete[] twiddle_lut;

	twiddle_lut_power = power;
	const uint count = 1 << power;
	twiddle_lut = new cfloat[count];

	const float arc = -2.0f*M_PI / (float)count;

	for (uint i = 0; i < count; i++){
		twiddle_lut[i].re = (fabsf(cos(arc*i)) > twiddle_epsilon)?(cos(arc*i)):(0.0f);
		twiddle_lut[i].im = (fabsf(sin(arc*i)) > twiddle_epsilon)?(sin(arc*i)):(0.0f);
	}
}

void freeTwiddleLUT(){
	delete[] twiddle_lut;
}

void fft(cfloat *out, cfloat *in, uint power){
	if (power == 0){
		*out = *in;
		return;
	}

	const uint count = 1 << (power-1);

	cfloat *even_in = new cfloat[count];
	cfloat *even_out = new cfloat[count];
	cfloat *odd_in = new cfloat[count];
	cfloat *odd_out = new cfloat[count];

	for (uint i = 0; i <  count; i++){
		even_in[i] = in[2*i];
		odd_in[i] = in[2*i + 1];
	}

	fft(even_out, even_in, power-1);
	fft(odd_out, odd_in, power-1);

	for (uint i = 0; i < count; i++){
		const cfloat t_even_out = 0.5f * even_out[i];
		const cfloat t_odd_out = 0.5f * odd_out[i] * twiddle_lut[i << (twiddle_lut_power-power)];
		out[i] = t_even_out + t_odd_out;
		out[i+count] = t_even_out - t_odd_out;
	}

	delete[] even_in;
	delete[] even_out;
	delete[] odd_in;
	delete[] odd_out;
}

void fft2(cfloat *out, cfloat *in, uint power){
	if (power == 0){
		*out = *in;
		return;
	}

	const uint dp = twiddle_lut_power - power;
	const uint ds = 1<<dp;
	const uint count = 1 << (power-1);

	fft2(out, in, power-1);
	fft2(out+count, in+ds, power-1);

	for (uint i = 0; i < count; i++){
		const cfloat tmp_l = 0.5f * out[i];
		const cfloat tmp_r = 0.5f * out[i + count] * twiddle_lut[i<<dp];

		out[i] = tmp_l + tmp_r;
		out[i+count] = tmp_l - tmp_r;
	}
}

void fft3(cfloat *out, cfloat *in, uint power){
	const uint back_shift = MAX_POWER - power;
	const uint count = 1 << power;

	// bit reverse ordering
	for (uint i = 0; i < count; i++)
		out[i] = in[bit_reverse(i) >> back_shift];

	// butterfly calculation
	for (uint m = 0; m < power; m++){
		const uint r_max = 1 << (power - m - 1);
		const uint shift = 1 << m;
		const uint dp = twiddle_lut_power - m -1;

		for(uint r = 0; r < r_max; r++){
			for (uint i = 0; i < shift; i++){
				const uint lidx = i + shift*(2*r);
				const uint ridx = i + shift*(2*r+1);

				const cfloat tmp_l = 0.5f * out[lidx];
				const cfloat tmp_r = 0.5f * out[ridx] * twiddle_lut[i<<dp];

				out[lidx] = tmp_l + tmp_r;
				out[ridx] = tmp_l - tmp_r;
			}
		}
	}
}

void fft4(cfloat *out, cfloat *in, uint power){
	const uint count = 1 << power;

	// bit reverse ordering
	for (uint i = 0; i < count; i++)
		out[i] = in[bit_reverse(i) >> (MAX_POWER - power)];

	// butterfly calculation
	for (uint m = 0; m < power; m++){
		const uint shift = 1 << m;
		const uint dp = twiddle_lut_power - m - 1;

		for (uint i = 0; i < count; ++i){
			if (shift & i)
				continue;

			const uint lidx = i;
			const uint ridx = i + shift;

			const cfloat tmp_l = 0.5f * out[lidx];
			const cfloat tmp_r = 0.5f * out[ridx] * twiddle_lut[(i%shift)<<dp];

			out[lidx] = tmp_l + tmp_r;
			out[ridx] = tmp_l - tmp_r;
		}
	}
}

int main(int argc, char const *argv[]){
	// cout << "sizeof(uint) = " << MAX_POWER << " bit" << endl;

	const uint power = 3;
	const uint count = 1 << power;

	initTwiddleLUT(power);

	cfloat *in = new cfloat[count];
	cfloat *out = new cfloat[count];
	
	for (uint i = 0; i < count; i++){
		in[i] = cfloat(i,0);
		show(in[i]);
	}
	cout << endl;

	fft4(out, in, power);

	for (uint i = 0; i < count; i++)
		show(out[i]);
	cout << endl;

	freeTwiddleLUT();

	uint x = 0x05;
	cout << x << ", ";
	x = x^0x04;
	cout << x << ", ";
	x = x^0x04;
	cout << x << endl;

	return 0;
}