#include "../microfacet.h"


Vector3 F_nu_metal(Real anisotropic, Real roughness, Real specularTint,
    Real specular, Real metallic, Real eta, const Frame& frame, const Vector3& baseColor, 
    const Vector3& dir_in, const Vector3& dir_out){

    Vector3 h = compute_h(dir_in, dir_out);
    Real aspect = sqrt(1.0f - 0.9f * anisotropic);

    Real alphaMin = 0.000001f;
    Real alphaX = std::max(alphaMin, pow(roughness, 2.0f) / aspect);
    Real alphaY = std::max(alphaMin, pow(roughness, 2.0f) * aspect);
    Vector3 hl = to_local(frame, h);

    Real DM = calculateD(alphaX, alphaY, hl);

    Vector3 dir_in_l = to_local(frame, dir_in);
    Vector3 dir_out_l = to_local(frame, dir_out);
    Real GM = calculateGM(dir_in_l, dir_out_l, alphaX, alphaY);

    Vector3 Ctint{1,1,1};
    Real lbc = luminance(baseColor);
    if(lbc > 0){
        Ctint = baseColor/lbc;
    }

    Vector3 Ks = (1-specularTint)*Vector3{1,1,1} + specularTint*Ctint;
    Vector3 C0 = specular*R0(eta)*(1-metallic)*Ks + metallic*baseColor;

    Real fac = 1-dot(h, dir_out);
    fac = pow(fac, 5.0);
    Vector3 Fm = C0 + (Vector3{1,1,1} - C0)*fac;

    Real constant = (.25/fabs(dot(frame.n, dir_in)));
    return constant*DM*GM*Fm;  
}

Real ImpSamp(PathVertex vertex, TransportDirection dir,
    Vector3 baseColor, Real specularTransmission, 
    Real metallic, Real subsurface, Real specular, Real roughness, 
    Real specularTint, Real anisotropic, Real sheen, Real sheenTint, 
    Real clearcoat, Real clearcoatGloss, Real eta, Frame frame, Vector3 dir_in, Vector3 dir_out){
    
    // Diffuse
    Real fdiff = (1-specularTransmission)*(1-metallic);
    Vector3 f_diff = F_diffuse(baseColor, frame.n, dir_in, dir_out, roughness, subsurface);
 
    // Sheen
    Real fsheen = (1-metallic)*sheen;
    Vector3 f_sheen = F_sheen(baseColor, sheenTint, frame.n, dir_in, dir_out);

    // Meta 
    Real fmet = (1-specularTransmission*(1-metallic));
    Vector3 f_metal = F_nu_metal(anisotropic, roughness, specularTint, specular,
        metallic, eta, frame, baseColor, dir_in, dir_out);

    // Clearcoat
    Real fcc = 0.25*clearcoat;
    Vector3 f_cc = F_cc(clearcoat, frame, dir_in, dir_out);
    
    if(dot(dir_in, vertex.geometric_normal) <= 0)
        f_metal = Vector3{0,0,0};

    // Glass
    Real fglass = (1-metallic)*specularTransmission;
    Vector3 f_glass = F_glass(anisotropic, roughness, 
        eta, vertex.geometric_normal, dir, frame,
        baseColor, dir_in, dir_out);

    if(dot(dir_in, vertex.geometric_normal) <= 0){
        f_diff = f_metal = f_cc = f_sheen = Vector3{0,0,0};
    }

    Vector3 f_disney = fdiff*f_diff + fsheen*f_sheen + fmet*f_metal + 
        fsheen*f_sheen + fglass*f_glass;

    Real sum = 0.333*(f_disney.x + f_disney.y + f_disney.z);

    return sum;

}

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;

    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Vector3 baseColor = eval(
    bsdf.base_color, 
    vertex.uv, 
    vertex.uv_screen_size, 
    texture_pool);

    Real specularTransmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specularTint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheenTint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = bsdf.eta;

    // Diffuse
    Real fdiff = (1-specularTransmission)*(1-metallic);
    Vector3 f_diff = F_diffuse(baseColor, frame.n, dir_in, dir_out, roughness, subsurface);
 
    // Sheen
    Real fsheen = (1-metallic)*sheen;
    Vector3 f_sheen = F_sheen(baseColor, sheenTint, frame.n, dir_in, dir_out);

    // Meta 
    Real fmet = (1-specularTransmission*(1-metallic));
    Vector3 f_metal = F_nu_metal(anisotropic, roughness, specularTint, specular,
        metallic, eta, frame, baseColor, dir_in, dir_out);

    // Clearcoat
    Real fcc = 0.25*clearcoat;
    Vector3 f_cc = F_cc(clearcoat, frame, dir_in, dir_out);
    
    if(dot(dir_in, vertex.geometric_normal) <= 0)
        f_metal = Vector3{0,0,0};

    // Glass
    Real fglass = (1-metallic)*specularTransmission;
    Vector3 f_glass = F_glass(anisotropic, roughness, 
        eta, vertex.geometric_normal, dir, frame,
        baseColor, dir_in, dir_out);

    if(dot(dir_in, vertex.geometric_normal) <= 0){
        f_diff = f_metal = f_cc = f_sheen = Vector3{0,0,0};
    }

    Vector3 f_disney = fdiff*f_diff + fsheen*f_sheen + fmet*f_metal + 
        fsheen*f_sheen + fglass*f_glass;

    return f_disney;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
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

    Vector3 baseColor = eval(
    bsdf.base_color, 
    vertex.uv, 
    vertex.uv_screen_size, 
    texture_pool);

    Real specularTransmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specularTint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheenTint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = bsdf.eta;

    Real f = ImpSamp(vertex, dir, baseColor, specularTransmission, metallic, subsurface, 
        specular, roughness, specularTint, anisotropic, sheen, sheenTint, clearcoat, 
        clearcoat_gloss, eta, frame, dir_in, dir_out);

    return f;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
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

    Vector3 baseColor = eval(
    bsdf.base_color, 
    vertex.uv, 
    vertex.uv_screen_size, 
    texture_pool);

    Real specularTransmission = eval(
        bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real metallic = eval(
        bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real subsurface = eval(
        bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specular = eval(
        bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real specularTint = eval(
        bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real anisotropic = eval(
        bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheen = eval(
        bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real sheenTint = eval(
        bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat = eval(
        bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real clearcoat_gloss = eval(
        bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real eta = bsdf.eta;

   // Diffuse
    Real fdiff = (1-specularTransmission)*(1-metallic);

    // Meta 
    Real fmet = (1-specularTransmission*(1-metallic));
 
    // Clearcoat
    Real fcc = 0.25*clearcoat;
 
    // Glass
    Real fglass = (1-metallic)*specularTransmission;

    std::vector<Real> wts = {fmet, fdiff, fglass, fcc};
    if(dot(vertex.geometric_normal, dir_in)) {
        wts = {0,0,1,0};
    }
    Real sum = 0;
    for(auto s : wts)
        sum += s;
    Real newsum = 0;
    for(auto& s : wts){
        s = s/sum;
        newsum += s;
    }

    std::vector<Real> bands = { wts[0], 
                                wts[0] + wts[1], 
                                wts[0] + wts[1] + wts[2],
                                1.0};

    if(dot(vertex.geometric_normal, dir_in) <= 0){

    }

    Real w = rnd_param_w; 
    if(w < bands[0]){
        // Metal
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
    } else if( w > bands[0] && w < bands[1]) {
        // Diffuse
               Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

        return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */};
    }
    else if(w > bands[1] && w < bands[2] ) {
        // GLASS
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

        Real w_new = (rnd_param_w - bands[1])  / (bands[2] - bands[1]);

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
    } else if (w > bands[2] && w < bands[3]) {
        // CLEARCOAT
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
}


TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
