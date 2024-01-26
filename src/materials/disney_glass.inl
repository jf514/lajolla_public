#include "../microfacet.h"

Real calcRS(Vector3 h, Real eta, Vector3 dir_in, Vector3 dir_out){
    Real num = dot(h, dir_in) - eta*dot(h, dir_out);
    Real den = dot(h, dir_in) + eta*dot(h, dir_out);

    return num/den;
}

Real calcRP(Vector3 h, Real eta, Vector3 dir_in, Vector3 dir_out){
    Real num = eta*dot(h, dir_in) - dot(h, dir_out);
    Real den = eta*dot(h, dir_in) + dot(h, dir_out);

    return num/den;
}

Real calcFg(Vector3 h, Real eta, Vector3 dir_in, Vector3 dir_out){
    Real RS = calcRS(h, eta, dir_in, dir_out);
    Real RP = calcRP(h, eta, dir_in, dir_out);

    return 0.5*(RS*RS + RP*RP);
}

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
 
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        h = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    Real Fg = calcFg(h, eta, dir_in, dir_out);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Vector3 baseColor = eval(
        bsdf.base_color, 
        vertex.uv, 
        vertex.uv_screen_size, 
        texture_pool);

    Real aspect = sqrt(1.0f - 0.9f * anisotropic);

    Real alphaMin = 0.000001f;
    Real alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    Real alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
    Vector3 hl = to_local(frame, h);
    Real DM = calculateD(alphaX, alphaY, hl);

    Vector3 dir_in_l = to_local(frame, dir_in);
    Vector3 dir_out_l = to_local(frame, dir_out);
    Real GM = calculateGM(dir_in_l, dir_out_l, alphaX, alphaY);

    Real h_dot_in = dot(h, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness) *
            smith_masking_gtr2(to_local(frame, dir_out), roughness);

    Vector3 res;
    if(reflect){
        res = F * DM * G * baseColor /(4 * fabs(dot(frame.n, dir_in)));
        //res = F * DM * G * baseColor /(4*fabs(dot(frame.n, dir_in)));
    } else {
        Vector3 baseColorSqrt(sqrt(baseColor.x), sqrt(baseColor.y), sqrt(baseColor.z));
        // Real num = (1-Fg) * DM * GM * fabs(dot(h, dir_out)*dot(h, dir_in));
        // Real den = fabs(dot(frame.n, dir_in))*pow((dot(h, dir_in) + eta * dot(h, dir_out)),2.);

        Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real h_dot_out = dot(h, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real fac =  (eta_factor * (1 - F) * DM * G * eta * eta * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
        res = fac * baseColorSqrt;
    }

    return res;
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    Real aspect = sqrt(1.0f - 0.9f * anisotropic);
    Real alphaMin = 0.000001f;
    Real alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    Real alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
    Vector3 hl = to_local(frame, half_vector);
    Real DM = calculateD(alphaX, alphaY, hl);


    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);



    Real D = GTR2(dot(half_vector, frame.n), roughness);
    Real G_in = smith_masking_gtr2(to_local(frame, dir_in), roughness);
    if (reflect) {
        return (F * DM * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * DM * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Now we need to decide whether to reflect or refract.
    // We do this using the Fresnel term.
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        // Reflection
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        // set eta to 0 since we are not transmitting
        return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
    } else {
        // Refraction
        // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        // (note that our eta is eta2 / eta1, and l = -dir_in)
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            return {};
        }
        // flip half_vector if needed
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out= sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{refracted, eta, roughness};
    }
}


// Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
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

// std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
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

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
