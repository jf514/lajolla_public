#include "../microfacet.h"


#include <cmath>

float R0(float eta) {
    return pow((eta - 1.0f) / (eta + 1.0f), 2.0f);
}

float calculateFc(Vector3 h, Vector3 dir_out) {
    Real R = R0(1.5);
    Real fac = pow(1.0f - std::abs(dot(h, dir_out)), 5.0f);
    Real Fc = R + (1-R)*fac;   
    return Fc;
}

float calculateDc(float alphaG, float hzl) {
    float algsq = alphaG*alphaG;
    float numerator = algsq - 1.0f;
    float denominator = c_PI * log(algsq) * (1.0f + numerator * hzl * hzl);

    return numerator / denominator;
}

// Function to calculate Lambda(omega)
float calculateLambdaCC(const Vector3& dir, float alphaX, float alphaY) {
    float numerator = pow((dir.x * alphaX), 2.0f) + pow((dir.y * alphaY), 2.0f);
    float denominator = pow(dir.z, 2.0f);
    return (sqrt(1.0f + numerator / denominator) - 1.0f) / 2.0f;
}

// Function to calculate G_m
float calculateGMCC(const Vector3& dir_in, const Vector3& dir_out, float alphaX, float alphaY) {

    float lambdaIn = calculateLambdaCC(dir_in, alphaX, alphaY);
    float lambdaOut = calculateLambdaCC(dir_out, alphaX, alphaY);

    float GIn = 1.0f / (1.0f + lambdaIn);
    float GOut = 1.0f / (1.0f + lambdaOut);

    return GIn * GOut;
}

Vector3 F_cc(Real ccg, const Frame& frame, const Vector3& dir_in, const Vector3& dir_out){
    Real alpha_g = (1-ccg)*0.1 + ccg*.001;

    Vector3 h = normalize(dir_in + dir_out);
    Real hzl = to_local(frame, h).z;
    Real Dc = calculateDc(alpha_g, hzl);

    Vector3 dir_in_l = to_local(frame, dir_in);
    Vector3 dir_out_l = to_local(frame, dir_out);
    Real Gc = calculateGMCC(dir_in_l, dir_out_l, 0.25, 0.25);
    Real Fc = calculateFc(h, dir_out);

    Vector3 res = (0.25/fabs(dot(frame.n, dir_in))) * Dc * Gc * Fc * Vector3(1,1,1);

    return res;
} 

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real ccg = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 res2 = F_cc(ccg, frame, dir_in, dir_out);

    return res2;
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return Real(0);
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

        // Homework 1: implement this!
    Real ccg = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real alpha_g = (1-ccg)*0.1 + ccg*.001;
    
    Vector3 h = normalize(dir_in + dir_out);
    Real hzl = to_local(frame, h).z;
    Real Dc = calculateDc(alpha_g, hzl);

    Real Gc = calculateGMCC(dir_in, dir_out, 0.25, 0.25);
    Real Fc = calculateFc(h, dir_out);

    Real res = (0.25/fabs(dot(frame.n, dir_in))) * Dc * fabs(dot(vertex.geometric_normal, h));

    // For Lambertian, we importance sample the cosine hemisphere domain.
    return res;
}

// Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
//     if (dot(vertex.geometric_normal, dir_in) < 0 ||
//             dot(vertex.geometric_normal, dir_out) < 0) {
//         // No light below the surface
//         return Real(0);
//     }
//     // Flip the shading frame if it is inconsistent with the geometry normal
//     Frame frame = vertex.shading_frame;
//     if (dot(frame.n, dir_in) < 0) {
//         frame = -frame;
//     }

//     // For Lambertian, we importance sample the cosine hemisphere domain.
//     return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
// }

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    // For Lambertian, we importance sample the cosine hemisphere domain.
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // Incoming direction is below the surface.
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

   Real ccg = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real alpha_g = (1-ccg)*0.1 + ccg*.001;
    Real alsq = alpha_g * alpha_g;

    Real u = rnd_param_uv[0];
    Real cos_hel = sqrt( (1 - pow(alsq, 1-u))/(1-alsq) );
    Real hel = acos(cos_hel);
    Real hazm = 2*c_PI*rnd_param_uv[1];
    Real hlx = sin(hel)*cos(hazm);
    Real hly = sin(hel)*sin(hazm);
    Real hlz = cos(hel);
    Vector3 hl{hlx, hly, hlz};
    hl = to_world(frame, hl);
    Vector3 reflected = normalize(-dir_in + 2. * dot(dir_in, hl)*hl);

    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, Real(1) /* roughness */};
}

// std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
//     // For Lambertian, we importance sample the cosine hemisphere domain.
//     if (dot(vertex.geometric_normal, dir_in) < 0) {
//         // Incoming direction is below the surface.
//         return {};
//     }
//     // Flip the shading frame if it is inconsistent with the geometry normal
//     Frame frame = vertex.shading_frame;
//     if (dot(frame.n, dir_in) < 0) {
//         frame = -frame;
//     }

//     return BSDFSampleRecord{
//         to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
//         Real(0) /* eta */, Real(1) /* roughness */};
// }

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
