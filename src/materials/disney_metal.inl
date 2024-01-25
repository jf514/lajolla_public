#include "../microfacet.h"


Vector3 compute_h(const Vector3& dir_in, const Vector3& dir_out){
    return normalize(dir_in+dir_out);
}

Vector3 calculateFM(const Vector3& baseColor, const Vector3& h, const Vector3& dir_out) {
    Real temp = 1.0f - fabs(dot(h, dir_out));
    temp = max(temp, Real(0));
    Real schlick =  pow(temp, 5.0);
    return baseColor + (Vector3{1.0, 1.0, 1.0} - baseColor) * schlick;
}

float calculateD(float alphaX, float alphaY, const Vector3& hl) {
    float hxl = hl.x / alphaX;
    float hyl = hl.y / alphaY;
    float hzl = hl.z;

    return 1.0f / (c_PI * alphaX * alphaY * pow(hxl * hxl + hyl * hyl + hzl * hzl, 2.0f));
}

// Function to calculate Lambda(omega)
float calculateLambda(const Vector3& dir, float alphaX, float alphaY) {
    float numerator = pow((dir.x * alphaX), 2.0f) + pow((dir.y * alphaY), 2.0f);
    float denominator = pow(dir.z, 2.0f);
    return (sqrt(1.0f + numerator / denominator) - 1.0f) / 2.0f;
}

// Function to calculate G_m
float calculateGM(const Vector3& dir_in, const Vector3& dir_out, float alphaX, float alphaY) {
    float lambdaIn = calculateLambda(dir_in, alphaX, alphaY);
    float lambdaOut = calculateLambda(dir_out, alphaX, alphaY);

    float GIn = 1.0f / (1.0f + lambdaIn);
    float GOut = 1.0f / (1.0f + lambdaOut);

    return GIn * GOut;
}

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Vector3 h = compute_h(dir_in, dir_out);

    Vector3 baseColor = eval(
        bsdf.base_color, 
        vertex.uv, 
        vertex.uv_screen_size, 
        texture_pool);

    Vector3 FM = calculateFM(baseColor, h, dir_out);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1.0f - 0.9f * anisotropic);

    Real alphaMin = 0.000001f;
    Real alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    Real alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
    Vector3 hl = to_local(frame, h);

    Real DM = calculateD(alphaX, alphaY, hl);

    Vector3 dir_in_l = to_local(frame, dir_in);
    Vector3 dir_out_l = to_local(frame, dir_out);
    Real GM = calculateGM(dir_in_l, dir_out_l, alphaX, alphaY);

    Vector3 res = (0.25/fabs(dot(frame.n, dir_in))) * DM * GM * FM;

    return res;
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {

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

    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    if (n_dot_out <= 0 || n_dot_h <= 0) {
        return 0;
    }

    // Spectrum S = eval(
    //     bsdf.specular_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Spectrum R = eval(
    //     bsdf.diffuse_reflectance, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Real lS = luminance(S), lR = luminance(R);
    // if (lS + lR <= 0) {
    //     return 0;
    // }
    Vector3 h = compute_h(dir_in, dir_out);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // We use the reflectance to determine whether to choose specular sampling lobe or diffuse.
    Real spec_prob = 1; //lS / (lS + lR);
    Real diff_prob = 0; //1 - spec_prob;
    // For the specular lobe, we use the ellipsoidal sampling from Heitz 2018
    // "Sampling the GGX Distribution of Visible Normals"
    // https://jcgt.org/published/0007/04/01/
    // this importance samples smith_masking(cos_theta_in) * GTR2(cos_theta_h, roughness) * cos_theta_out
    //Real G = calculateGM(dir_in, dir_out, alphaX, alphaY);
    //Real D = GTR2(n_dot_h, roughness);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1.0f - 0.9f * anisotropic);


    Real alphaMin = 0.000001f;
    Real alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    Real alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
    Vector3 hl = to_local(frame, h);

    Real DM = calculateD(alphaX, alphaY, hl);

    Vector3 dir_in_l = to_local(frame, dir_in);
    Vector3 dir_out_l = to_local(frame, dir_out);
    Real GM = calculateGM(dir_in_l, dir_out_l, alphaX, alphaY);


    // (4 * cos_theta_v) is the Jacobian of the reflectiokn
    spec_prob *= (GM * DM) / (4 * n_dot_in);
    // For the diffuse lobe, we importance sample cos_theta_out
    diff_prob *= n_dot_out / c_PI;
    return spec_prob + diff_prob;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // Sample from the specular lobe.
    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
    
    // Transform the micro normal to world space
    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Reflect over the world space normal
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        Real(0) /* eta */, roughness /* roughness */
    };
 
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
