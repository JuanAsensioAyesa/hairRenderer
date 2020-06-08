
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// tests/hair.cpp*
#include <yocto/yocto_math.h>
#include <yocto_extension/spectrum.h>
// #include <yocto_extension/yocto_extension.h>

#include <iostream>

namespace math = yocto::math;
namespace ex   = yocto::extension;
using math::rand1f;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;
// #include "materials/hair.h"
// #include "rng.h"
// #include "sampling.h"
// #include "tests/gtest/gtest.h"

// using namespace pbrt;

// Hair Tests
vec3f UniformSampleSphere(const vec2f &u) {
  float z   = 1 - 2 * u[0];
  float r   = std::sqrt(std::max((float)0, (float)1 - z * z));
  float phi = 2 * math::pi * u[1];
  return vec3f(r * std::cos(phi), r * std::sin(phi), z);
}

static float Inv4Pi2 = 0.07957747154594766788;
float        UniformSpherePdf() { return Inv4Pi2; }

// void tryReciprocity() {
//   yocto::math::rng_state rng = math::rng_state();
//   //   float random = rand1f(rng);
//   for (int i = 0; i < 10; ++i) {
//     Hair h(-1 + 2 * rng.Uniformfloat(), 1.55,
//            HairBSDF::SigmaAFromConcentration(.3 + 7.7 * rng.Uniformfloat()),
//            .1 + .9 * rng.Uniformfloat(),
//            .1 + .9 * rng.Uniformfloat());

//     vec3f          wi = UniformSampleSphere({rand1f(rng), rand1f(rng)});
//     vec3f          wo = UniformSampleSphere({rand1f(rng), rand1f(rng)});
//     pbrt::Spectrum a  = yocto::extension::f(wi, wo) *
//                        yocto::extension::AbsCosTheta(wo);
//     pbrt::Spectrum b = yocto::extension::f(wo, wi) *
//                        yocto::extension::AbsCosTheta(wi);

//     if (a == b) {
//       std::cout << "Son iguales" << std::endl;
//     } else {
//       std::cout << "NO" << std::endl;
//       std::cout << a << std::endl;
//       std::cout << b << std::endl;
//       std::cout << wi.x << " " << wi.y << " " << wi.z << std::endl;
//       std::cout << wo.x << " " << wo.y << " " << wo.z << std::endl;
//     }
//     // EXPECT_EQ(a.y(), b.y()) << h << ", a = " << a << ", b = " << b
//     //                         << ", wi = " << wi << ", wo = " << wo;
//   }
//   int n;
//   std::cin >> n;
// }

// TEST(Hair, WhiteFurnace) {
//   RNG      rng;
//   Vector3f wo = UniformSampleSphere({rng.Uniformfloat(),
//   rng.Uniformfloat()}); for (float beta_m = .1; beta_m < 1; beta_m += .2) {
//     for (float beta_n = .1; beta_n < 1; beta_n += .2) {
//       // Estimate reflected uniform incident radiance from hair
//       Spectrum sum   = 0.f;
//       int      count = 300000;
//       for (int i = 0; i < count; ++i) {
//         float    h       = -1 + 2. * rng.Uniformfloat();
//         Spectrum sigma_a = 0.f;
//         HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);
//         Vector3f wi = UniformSampleSphere(
//             {rng.Uniformfloat(), rng.Uniformfloat()});
//         sum += hair.f(wo, wi) * AbsCosTheta(wi);
//       }
//       float avg = sum.y() / (count * UniformSpherePdf());
//       EXPECT_TRUE(avg >= .95 && avg <= 1.05);
//     }
//   }
// }

// TEST(Hair, WhiteFurnaceSampled) {
//   RNG      rng;
//   Vector3f wo = UniformSampleSphere({rng.Uniformfloat(),
//   rng.Uniformfloat()}); for (float beta_m = .1; beta_m < 1; beta_m += .2) {
//     for (float beta_n = .1; beta_n < 1; beta_n += .2) {
//       Spectrum sum   = 0.f;
//       int      count = 300000;
//       for (int i = 0; i < count; ++i) {
//         float    h       = -1 + 2. * rng.Uniformfloat();
//         Spectrum sigma_a = 0.f;
//         HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);

//         Vector3f wi;
//         float    pdf;
//         Point2f  u = {rng.Uniformfloat(), rng.Uniformfloat()};
//         Spectrum f = hair.Sample_f(wo, &wi, u, &pdf, nullptr);
//         if (pdf > 0) sum += f * AbsCosTheta(wi) / pdf;
//       }
//       float avg = sum.y() / count;
//       EXPECT_TRUE(avg >= .99 && avg <= 1.01) << avg;
//     }
//   }
// }

// TEST(Hair, SamplingWeights) {
//   RNG rng;
//   for (float beta_m = .1; beta_m < 1; beta_m += .2)
//     for (float beta_n = .4; beta_n < 1; beta_n += .2) {
//       int count = 10000;
//       for (int i = 0; i < count; ++i) {
//         // Check _HairBSDF::Sample\_f()_ sample weight
//         float    h       = -1 + 2 * rng.Uniformfloat();
//         Spectrum sigma_a = 0;
//         HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);
//         Vector3f wo = UniformSampleSphere(
//             {rng.Uniformfloat(), rng.Uniformfloat()});
//         Vector3f wi;
//         float    pdf;
//         Point2f  u = {rng.Uniformfloat(), rng.Uniformfloat()};
//         Spectrum f = hair.Sample_f(wo, &wi, u, &pdf, nullptr);
//         if (pdf > 0) {
//           // Verify that hair BSDF sample weight is close to 1 for
//           // _wi_
//           EXPECT_GT(f.y() * AbsCosTheta(wi) / pdf, 0.999);
//           EXPECT_LT(f.y() * AbsCosTheta(wi) / pdf, 1.001);
//         }
//       }
//     }
// }

// TEST(Hair, SamplingConsistency) {
//   RNG rng;
//   for (float beta_m = .2; beta_m < 1; beta_m += .2)
//     for (float beta_n = .4; beta_n < 1; beta_n += .2) {
//       // Declare variables for hair sampling test
//       const int count   = 64 * 1024;
//       Spectrum  sigma_a = .25;
//       Vector3f  wo      = UniformSampleSphere(
//           {rng.Uniformfloat(), rng.Uniformfloat()});
//       auto     Li = [](const Vector3f &w) -> Spectrum { return w.z * w.z; };
//       Spectrum fImportance = 0, fUniform = 0;
//       for (int i = 0; i < count; ++i) {
//         // Compute estimates of scattered radiance for hair sampling
//         // test
//         float    h = -1 + 2 * rng.Uniformfloat();
//         HairBSDF hair(h, 1.55, sigma_a, beta_m, beta_n, 0.f);
//         Vector3f wi;
//         float    pdf;
//         Point2f  u = {rng.Uniformfloat(), rng.Uniformfloat()};
//         Spectrum f = hair.Sample_f(wo, &wi, u, &pdf, nullptr);
//         if (pdf > 0)
//           fImportance += f * Li(wi) * AbsCosTheta(wi) / (count * pdf);
//         wi = UniformSampleSphere(u);
//         fUniform += hair.f(wo, wi) * Li(wi) * AbsCosTheta(wi) /
//                     (count * UniformSpherePdf());
//       }
//       // Verify consistency of estimated hair reflected radiance values
//       float err = std::abs(fImportance.y() - fUniform.y()) / fUniform.y();
//       EXPECT_LT(err, 0.05);
//     }
// }
