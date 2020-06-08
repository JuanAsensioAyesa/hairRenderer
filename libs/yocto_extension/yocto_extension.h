//
// # Yocto/Extension: Tiny Yocto/GL extension
//
//

//
// LICENSE:
//
// Copyright (c) 2020 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_EXTENSION_H_
#define _YOCTO_EXTENSION_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

#include <atomic>
#include <future>
#include <memory>
#include <numeric>

#include "pbrt.h"
#include "spectrum.h"
#include "stringprint.h"

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::extension {

// Namespace aliases
namespace ext = yocto::extension;
namespace img = yocto::image;

// Math defitions
using math::abs;
using math::acos;
using math::atan2;
using math::bbox3f;
using math::byte;
using math::clamp;
using math::cos;
using math::exp;
using math::flt_max;
using math::fmod;
using math::frame3f;
using math::fresnel_conductor;
using math::fresnel_dielectric;
using math::identity3x3f;
using math::identity3x4f;
using math::invalidb3f;
using math::log;
using math::make_rng;
using math::min;
using math::pif;
using math::pow;
using math::radians;
using math::ray3f;
using math::rng_state;
using math::sample_discrete_cdf;
using math::sample_discrete_cdf_pdf;
using math::sample_uniform;
using math::sample_uniform_pdf;
using math::sin;
using math::sqrt;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;
using math::zero2f;
using math::zero2i;
using math::zero3f;
using math::zero3i;
using math::zero4f;
using math::zero4i;

}  // namespace yocto::extension

// -----------------------------------------------------------------------------
// HIGH LEVEL API
// -----------------------------------------------------------------------------
namespace yocto::extension {

// HairBSDF Constants
static const int   pMax        = 3;
static const float SqrtPiOver8 = 0.626657069f;

// HairBSDF Declarations
class HairBSDF {
 public:
  // HairBSDF Public Methods
  HairBSDF(float h, float eta, const pbrt::Spectrum &sigma_a, float beta_m,
      float beta_n, float alpha);
  pbrt::Spectrum f(const vec3f &wo, const vec3f &wi) const;
  pbrt::Spectrum Sample_f(
      const vec3f &wo, vec3f *wi, const vec2f &u, float *pdf) const;
  float                 Pdf(const vec3f &wo, const vec3f &wi) const;
  std::string           ToString() const;
  static pbrt::Spectrum SigmaAFromConcentration(float ce, float cp);
  static pbrt::Spectrum SigmaAFromReflectance(
      const pbrt::Spectrum &c, float beta_n);

 private:
  // HairBSDF Private Methods
  std::array<float, pMax + 1> ComputeApPdf(float cosThetaO) const;

  // HairBSDF Private Data
  float          h, gammaO, eta;
  pbrt::Spectrum sigma_a;
  float          beta_m, beta_n;
  float          v[pMax + 1];
  float          s;
  float          sin2kAlpha[3], cos2kAlpha[3];
};

// Refraction index of the interior of the hair
// float eta = 1.55;

// number of transmissions

// float sigmaAux[3] = {0.06, 0.1, 0.2};

// // Absortion coefficient of the interior of the hair
// // vec3f sigma_a = vec3f{0.06, 0.1, 0.2};
// pbrt::Spectrum sigma_a = pbrt::Spectrum::FromRGB(sigmaAux);
// // pbrt::Spectrum sigma_a = 0.25;
// float alpha  = 2.0;
// float h      = 0.43; /* TODO: hay que cambiar el 0.5 por algo nandom*/
// float beta_n = 0.3;
// float beta_m = 0.125;

inline float Sqr(float v) { return v * v; }
// efficient Pow function
template <int n>
static float Pow(float v) {
  static_assert(n > 0, "Power can’t be negative");
  float n2 = Pow<n / 2>(v);
  return n2 * n2 * Pow<n & 1>(v);
}

template <>
float Pow<1>(float v) {
  return v;
}

template <>
float Pow<0>(float v) {
  return 1;
}

//
inline float SafeASin(float x) {
  if (x >= -1.0001 && x <= 1.0001)
    return std::asin(yocto::math::clamp(x, -1.0, 1.0));
}

//
inline float SafeSqrt(float x) {
  if (x >= -1e-4) {
    return std::sqrt(std::max(float(0), x));
  }
}
inline float AbsCosTheta(const vec3f &w) { return std::abs(w.z); }

inline float I0(float x), LogI0(float x);

// Hair Local Functions
static float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
    float sinThetaO, float v) {
  float a = cosThetaI * cosThetaO / v;
  float b = sinThetaI * sinThetaO / v;
  float mp =
      (v <= .1)
          ? (std::exp(LogI0(a) - b - 1 / v + 0.6931f + std::log(1 / (2 * v))))
          : (std::exp(-b) * I0(a)) / (std::sinh(1 / v) * 2 * v);
  // CHECK(!std::isinf(mp) && !std::isnan(mp));

  return mp;
}

inline float I0(float x) {
  float   val   = 0;
  float   x2i   = 1;
  int64_t ifact = 1;
  int     i4    = 1;
  // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
  for (int i = 0; i < 10; ++i) {
    if (i > 1) ifact *= i;
    val += x2i / (i4 * Sqr(ifact));
    x2i *= x * x;
    i4 *= 4;
  }
  return val;
}

inline float LogI0(float x) {
  if (x > 12)
    return x + 0.5 * (-std::log(2 * yocto::math::pi) + std::log(1 / x) +
                         1 / (8 * x));
  else
    return std::log(I0(x));
}

inline float Phi(int p, float gammaO, float gammaT) {
  return 2 * p * gammaT - 2 * gammaO + p * yocto::math::pi;
}

inline float Logistic(float x, float s) {
  x = std::abs(x);
  return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}

inline float LogisticCDF(float x, float s) {
  return 1 / (1 + std::exp(-x / s));
}

inline float TrimmedLogistic(float x, float s, float a, float b) {
  // CHECK_LT(a, b);
  return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

static float SampleTrimmedLogistic(float u, float s, float a, float b) {
  // CHECK_LT(a, b);
  float k = LogisticCDF(b, s) - LogisticCDF(a, s);
  float x = -s * std::log(1 / (u * k + LogisticCDF(a, s)) - 1);
  // CHECK(!std::isnan(x));
  return clamp(x, a, b);
}

inline float Np(float phi, int p, float s, float gammaO, float gammaT) {
  float dphi = phi - Phi(p, gammaO, gammaT);
  // Remap _dphi_ to $[-\pi,\pi]$
  while (dphi > yocto::math::pi) dphi -= 2 * yocto::math::pi;
  while (dphi < -yocto::math::pi) dphi += 2 * yocto::math::pi;
  return TrimmedLogistic(dphi, s, -yocto::math::pi, yocto::math::pi);
}

float FrDielectric(float cosThetaI, float etaI, float etaT) {
  cosThetaI = clamp(cosThetaI, -1.0, 1.0);
  // Potentially swap indices of refraction
  bool entering = cosThetaI > 0.f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = std::abs(cosThetaI);
  }

  // Compute _cosThetaT_ using Snell's law
  float sinThetaI = std::sqrt(std::max((float)0, 1 - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;

  // Handle total internal reflection
  if (sinThetaT >= 1) return 1;
  float cosThetaT = std::sqrt(std::max((float)0, 1 - sinThetaT * sinThetaT));
  float Rparl     = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                ((etaT * cosThetaI) + (etaI * cosThetaT));
  float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                ((etaI * cosThetaI) + (etaT * cosThetaT));
  return (Rparl * Rparl + Rperp * Rperp) / 2;
}

// static std::array<pbrt::Spectrum, pMax + 1> Ap(
//     float cosThetaO, float eta, float h, const pbrt::Spectrum &T) {
//   std::array<pbrt::Spectrum, pMax + 1> ap;
//   // Compute $p=0$ attenuation at initial cylinder intersection
//   float cosGammaO = SafeSqrt(1 - h * h);
//   float cosTheta  = cosThetaO * cosGammaO;
//   float f         = FrDielectric(cosTheta, 1.f, eta);

//   // std::cout << "Dielectric " << f << std::endl;
//   ap[0] = f;

//   // Compute $p=1$ attenuation term
//   ap[1] = sqrt(1 - f) * T;

//   // Compute attenuation terms up to $p=_pMax_$
//   for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

//   // Compute attenuation term accounting for remaining orders of scattering
//   ap[pMax] = ap[pMax - 1] * f * T / (pbrt::Spectrum(1.f) - T * f);
//   return ap;
// }

// HairBSDF Method Definitions
HairBSDF::HairBSDF(float h2, float eta2, const pbrt::Spectrum &sigma_a2,
    float beta_m2, float beta_n2, float alpha2) {
  h       = h2;
  gammaO  = SafeASin(h2);
  eta     = eta2;
  sigma_a = sigma_a2;
  beta_m  = beta_m2;
  beta_n  = beta_n2;
  // CHECK(h >= -1 && h <= 1);
  // CHECK(beta_m >= 0 && beta_m <= 1);
  // CHECK(beta_n >= 0 && beta_n <= 1);
  // Compute longitudinal variance from $\beta_m$
  // static_assert(pMax >= 3,
  //     "Longitudinal variance code must be updated to handle low pMax");
  v[0] = Sqr(0.726f * beta_m + 0.812f * Sqr(beta_m) + 3.7f * Pow<20>(beta_m));
  v[1] = .25 * v[0];
  v[2] = 4 * v[0];
  for (int p = 3; p <= pMax; ++p)
    // TODO: is there anything better here?
    v[p] = v[2];

  // Compute azimuthal logistic scale factor from $\beta_n$
  s = SqrtPiOver8 *
      (0.265f * beta_n + 1.194f * Sqr(beta_n) + 5.372f * Pow<22>(beta_n));
  // CHECK(!std::isnan(s));

  // Compute $\alpha$ terms for hair scales
  sin2kAlpha[0] = std::sin(math::radians(alpha2));
  cos2kAlpha[0] = SafeSqrt(1 - Sqr(sin2kAlpha[0]));
  for (int i = 1; i < 3; ++i) {
    sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
    cos2kAlpha[i] = Sqr(cos2kAlpha[i - 1]) - Sqr(sin2kAlpha[i - 1]);
  }
}

static uint32_t Compact1By1(uint32_t x) {
  // TODO: as of Haswell, the PEXT instruction could do all this in a
  // single instruction.
  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x &= 0x55555555;
  // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 1)) & 0x33333333;
  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;
  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;
  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;
  return x;
}

static vec2f Demuxfloat(float f) {
  // CHECK(f >= 0 && f < 1);
  uint64_t v = f * (1ull << 32);
  // CHECK_LT(v, 0x100000000);
  uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
  return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
}

static std::array<pbrt::Spectrum, pMax + 1> Ap(
    float cosThetaO, float eta, float h, const pbrt::Spectrum &T) {
  std::array<pbrt::Spectrum, pMax + 1> ap;
  // Compute $p=0$ attenuation at initial cylinder intersection
  float cosGammaO = SafeSqrt(1 - h * h);
  float cosTheta  = cosThetaO * cosGammaO;
  float f         = FrDielectric(cosTheta, 1.f, eta);
  ap[0]           = f;

  // Compute $p=1$ attenuation term
  ap[1] = Sqr(1 - f) * T;

  // Compute attenuation terms up to $p=_pMax_$
  for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;

  // Compute attenuation term accounting for remaining orders of scattering
  ap[pMax] = ap[pMax - 1] * f * T / (pbrt::Spectrum(1.f) - T * f);
  return ap;
}

pbrt::Spectrum HairBSDF::f(const vec3f &wo, const vec3f &wi) const {
  // Compute hair coordinate system terms related to _wo_
  float sinThetaO = wo.x;
  float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  float phiO      = std::atan2(wo.z, wo.y);

  // Compute hair coordinate system terms related to _wi_
  float sinThetaI = wi.x;
  float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  float phiI      = std::atan2(wi.z, wi.y);

  // Compute $\cos \thetat$ for refracted ray
  float sinThetaT = sinThetaO / eta;
  float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

  // Compute $\gammat$ for refracted ray
  float etap      = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
  float gammaT    = SafeASin(sinGammaT);

  // Compute the transmittance _T_ of a single path through the cylinder
  pbrt::Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));

  // Evaluate hair BSDF
  float                                phi = phiI - phiO;
  std::array<pbrt::Spectrum, pMax + 1> ap  = Ap(cosThetaO, eta, h, T);
  pbrt::Spectrum                       fsum(0.);
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = std::abs(cosThetaOp);
    fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * ap[p] *
            Np(phi, p, s, gammaO, gammaT);
  }

  // Compute contribution of remaining terms after _pMax_
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * ap[pMax] /
          (2.f * yocto::math::pi);
  if (AbsCosTheta(wi) > 0) fsum /= AbsCosTheta(wi);
  // CHECK(!std::isinf(fsum.y()) && !std::isnan(fsum.y()));
  return fsum;
}

std::array<float, pMax + 1> HairBSDF::ComputeApPdf(float cosThetaO) const {
  // Compute array of $A_p$ values for _cosThetaO_
  float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);

  // Compute $\cos \thetat$ for refracted ray
  float sinThetaT = sinThetaO / eta;
  float cosThetaT = SafeSqrt(1 - Sqr(sinThetaT));

  // Compute $\gammat$ for refracted ray
  float etap      = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));

  // Compute the transmittance _T_ of a single path through the cylinder
  pbrt::Spectrum T = Exp(-sigma_a * (2 * cosGammaT / cosThetaT));
  std::array<pbrt::Spectrum, pMax + 1> ap = Ap(cosThetaO, eta, h, T);

  // Compute $A_p$ PDF from individual $A_p$ terms
  std::array<float, pMax + 1> apPdf;
  float sumY = std::accumulate(ap.begin(), ap.end(), float(0),
      [](float s, const pbrt::Spectrum &ap) { return s + ap.y(); });
  for (int i = 0; i <= pMax; ++i) apPdf[i] = ap[i].y() / sumY;
  return apPdf;
}

pbrt::Spectrum HairBSDF::Sample_f(
    const vec3f &wo, vec3f *wi, const vec2f &u2, float *pdf) const {
  // Compute hair coordinate system terms related to _wo_
  float sinThetaO = wo.x;
  float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  float phiO      = std::atan2(wo.z, wo.y);

  // Derive four random samples from _u2_
  vec2f u[2] = {Demuxfloat(u2[0]), Demuxfloat(u2[1])};

  // Determine which term $p$ to sample for hair scattering
  std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);
  int                         p;
  for (p = 0; p < pMax; ++p) {
    if (u[0][0] < apPdf[p]) break;
    u[0][0] -= apPdf[p];
  }

  // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
  float sinThetaOp, cosThetaOp;
  if (p == 0) {
    sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
    cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
  } else if (p == 1) {
    sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
    cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
  } else if (p == 2) {
    sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
    cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
  } else {
    sinThetaOp = sinThetaO;
    cosThetaOp = cosThetaO;
  }

  // Sample $M_p$ to compute $\thetai$
  u[1][0] = std::max(u[1][0], float(1e-5));
  float cosTheta =
      1 + v[p] * std::log(u[1][0] + (1 - u[1][0]) * std::exp(-2 / v[p]));
  float sinTheta  = SafeSqrt(1 - Sqr(cosTheta));
  float cosPhi    = std::cos(2 * yocto::math::pi * u[1][1]);
  float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));

  // Sample $N_p$ to compute $\Delta\phi$

  // Compute $\gammat$ for refracted ray
  float etap      = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float gammaT    = SafeASin(sinGammaT);
  float dphi;
  if (p < pMax)
    dphi = Phi(p, gammaO, gammaT) +
           SampleTrimmedLogistic(u[0][1], s, -yocto::math::pi, yocto::math::pi);
  else
    dphi = 2 * yocto::math::pi * u[0][1];

  // Compute _wi_ from sampled hair scattering angles
  float phiI = phiO + dphi;
  *wi        = vec3f(
      sinThetaI, cosThetaI * std::cos(phiI), cosThetaI * std::sin(phiI));

  // Compute PDF for sampled hair scattering direction _wi_
  *pdf = 0;
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = std::abs(cosThetaOp);
    *pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * apPdf[p] *
            Np(dphi, p, s, gammaO, gammaT);
  }
  *pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) *
          apPdf[pMax] * (1 / (2 * yocto::math::pi));
  // if (std::abs(wi->x) < .9999) CHECK_NEAR(*pdf, Pdf(wo, *wi), .01);
  return f(wo, *wi);
}

float HairBSDF::Pdf(const vec3f &wo, const vec3f &wi) const {
  // Compute hair coordinate system terms related to _wo_
  float sinThetaO = wo.x;
  float cosThetaO = SafeSqrt(1 - Sqr(sinThetaO));
  float phiO      = std::atan2(wo.z, wo.y);

  // Compute hair coordinate system terms related to _wi_
  float sinThetaI = wi.x;
  float cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  float phiI      = std::atan2(wi.z, wi.y);

  // Compute $\gammat$ for refracted ray
  float etap      = std::sqrt(eta * eta - Sqr(sinThetaO)) / cosThetaO;
  float sinGammaT = h / etap;
  float gammaT    = SafeASin(sinGammaT);

  // Compute PDF for $A_p$ terms
  std::array<float, pMax + 1> apPdf = ComputeApPdf(cosThetaO);

  // Compute PDF sum for hair scattering events
  float phi = phiI - phiO;
  float pdf = 0;
  for (int p = 0; p < pMax; ++p) {
    // Compute $\sin \thetao$ and $\cos \thetao$ terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * cos2kAlpha[1] - cosThetaO * sin2kAlpha[1];
      cosThetaOp = cosThetaO * cos2kAlpha[1] + sinThetaO * sin2kAlpha[1];
    }

    // Handle remainder of $p$ values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * cos2kAlpha[0] + cosThetaO * sin2kAlpha[0];
      cosThetaOp = cosThetaO * cos2kAlpha[0] - sinThetaO * sin2kAlpha[0];
    } else if (p == 2) {
      sinThetaOp = sinThetaO * cos2kAlpha[2] + cosThetaO * sin2kAlpha[2];
      cosThetaOp = cosThetaO * cos2kAlpha[2] - sinThetaO * sin2kAlpha[2];
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }

    // Handle out-of-range $\cos \thetao$ from scale adjustment
    cosThetaOp = std::abs(cosThetaOp);
    pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, v[p]) * apPdf[p] *
           Np(phi, p, s, gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, v[pMax]) * apPdf[pMax] *
         (1 / (2 * yocto::math::pi));
  return pdf;
}

// std::string HairBSDF::ToString() const {
//   return StringPrintf(
//              "[ Hair h: %f gammaO: %f eta: %f beta_m: %f beta_n: %f "
//              "v[0]: %f s: %f sigma_a: ",
//              h, gammaO, eta, beta_m, beta_n, v[0], s) +
//          sigma_a.ToString() + std::string("  ]");
// }

pbrt::Spectrum HairBSDF::SigmaAFromConcentration(float ce, float cp) {
  float sigma_a[3];
  float eumelaninSigmaA[3]   = {0.419f, 0.697f, 1.37f};
  float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
  for (int i = 0; i < 3; ++i)
    sigma_a[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
  return pbrt::Spectrum::FromRGB(sigma_a);
}

pbrt::Spectrum HairBSDF::SigmaAFromReflectance(
    const pbrt::Spectrum &c, float beta_n) {
  pbrt::Spectrum sigma_a;
  for (int i = 0; i < pbrt::Spectrum::nSamples; ++i)
    sigma_a[i] = Sqr(
        std::log(c[i]) / (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
                             10.73f * Pow<3>(beta_n) + 5.574f * Pow<4>(beta_n) +
                             0.245f * Pow<5>(beta_n)));
  return sigma_a;
}

}  // namespace yocto::extension

//#include "tests.h"

#endif
