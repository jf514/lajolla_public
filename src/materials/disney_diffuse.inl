Real FD_90(const Vector3 &dir_out, const Vector3& h, Real roughness){
    Real F_D90 = (0.5 + 2 * roughness * pow(dot(h, dir_out),2));
    return F_D90;
} 

Real FD(const Vector3 &norm, 
    const Vector3 &dir, 
    const Vector3 &dir_out, 
    const Vector3& h, 
    Real roughness){
    Real FD90 = FD_90(dir_out, h, roughness);
    Real FD = (1. + (FD90 - 1.)*pow(1. - fabs(dot(norm, dir)),5));
    return FD;
}

Real FSS90(const Vector3& dir_out, const Vector3& h, Real roughness){
    return roughness * pow(fabs(dot(h, dir_out)),2);
}

Real FSS(const Vector3& norm, const Vector3& dir, 
    const Vector3& dir_out, const Vector3& h, Real roughness){
        Real FSS_90 = FSS90(dir_out, h, roughness);
        Real FSS = (1 + 
                        (FSS_90 - 1) * 
                        (pow( max(1 - dot(norm, dir),Real(0)), 5))
                    );
        return FSS;
}

Real FSS_Combined(const Vector3& norm, const Vector3& dir_in, const Vector3& dir_out,
    const Vector3& h, Real roughness){
        Real FSS_in = FSS(norm, dir_in, dir_out, h, roughness);
        Real FSS_out = FSS(norm, dir_out, dir_out, h, roughness);
        Real omega_term = (1./(fabs(dot(norm, dir_in)) + fabs(dot(norm, dir_out)))); 
        Real fac = (FSS_in * FSS_out*(omega_term - 0.5)+0.5)*fabs(dot(norm,dir_out));
        return fac;
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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

    // Base diffuse
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        
    Vector3 h = normalize(dir_in + dir_out);
    Real FD_in = FD(frame.n, dir_in, dir_out, h, roughness);
    Real FD_out = FD(frame.n, dir_out, dir_out, h, roughness);
    
    Real f_base_diffuse = (FD_in * FD_out/c_PI) * 
        fabs(dot(frame.n, dir_out));

    // Subsurface
    Real FSS_comb = FSS_Combined(frame.n, dir_in, dir_out, h, roughness);
    Real f_sub_surf = (1.25/c_PI)*FSS_comb; 
   
    Real subsurface = eval(
        bsdf.subsurface, 
        vertex.uv, 
        vertex.uv_screen_size,
        texture_pool
    );

    Vector3 baseColor = eval(
        bsdf.base_color, 
        vertex.uv, 
        vertex.uv_screen_size, 
        texture_pool);

    Real diffuse = (1-subsurface)*f_base_diffuse + subsurface*f_sub_surf;

    return diffuse * baseColor;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // Importance sample the hemisphere
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
        Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, roughness /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
