#include "../microfacet.h"


Spectrum eval_op::operator()(const DisneySheen &bsdf) const {
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

    Vector3 baseColor = eval(
        bsdf.base_color, 
        vertex.uv, 
        vertex.uv_screen_size, 
        texture_pool);

    Real sheenTint = eval(
        bsdf.sheen_tint, 
        vertex.uv, 
        vertex.uv_screen_size, 
        texture_pool);

    Vector3 Ctint{1,1,1};
    Real lbc = luminance(baseColor);
    if(lbc > 0){
        Ctint = baseColor/lbc;
    }

    Vector3 Csheen = (1-sheenTint)*Vector3(1,1,1) + sheenTint*Ctint;
    Vector3 h = normalize(dir_in + dir_out);
    Real fac = pow(1.0f - std::abs(dot(h, dir_out)), 5.0f)*fabs(dot(frame.n, dir_out));
    Vector3 fsheen = fac*Csheen;

    // Homework 1: implement this!
    return fsheen;
}

Real pdf_sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
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

    // For DisneySheen, we importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneySheen &bsdf) const {
    // For DisneySheen, we importance sample the cosine hemisphere domain.
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // Incoming direction is below the surface.
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}


TextureSpectrum get_texture_op::operator()(const DisneySheen &bsdf) const {
    return bsdf.base_color;
}
