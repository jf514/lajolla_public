#include "../microfacet.h"


#include <cmath>


float R0(float eta) {
    return pow((eta - 1.0f) / (eta + 1.0f), 2.0f);
}

float calculateFc(Vector3 h, Vector3 dir_out) {
    Real R = R0(1.5);
    Real fac = pow(1.0f - std::abs(dot(h, dir_out)), 5.0f);
    Real Fc = R + (1-R*fac);   
    return Fc;
}

float calculateDc(float alphaG, float hzl) {
    float numerator = pow(alphaG, 2.0f) - 1.0f;
    float denominator = c_PI * log(pow(alphaG, 2.0f)) * (1.0f + numerator * pow(hzl, 2.0f));

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

    Real alpha_g = (1-ccg)*0.1 + ccg*.001;

    Vector3 h = normalize(dir_in + dir_out);
    Real hzl = to_local(frame, h).z;
    Real Dc = calculateDc(alpha_g, hzl);

    Real Gc = calculateGMCC(dir_in, dir_out, 0.25, 0.25);
    Real Fc = calculateFc(h, dir_out);

    Vector3 res = (0.25/fabs(dot(frame.n, dir_in))) * Dc * Gc * Fc * Vector3(1,1,1);

    return res;
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    return 0;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
 


    return {};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
