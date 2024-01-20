Real FD_90(const Vector3 &dir_out, const Vector3& h, Real roughness){
    Real F_D90 = (0.5 + 2 * roughness * dot(h, dir_out)*dot(h, dir_out));
    return F_D90;
} 

Real FD(const Vector3 &norm, const Vector3 &dir, const Vector3 &dir_out, const Vector3& h, Real roughness){
    Real FD90 = FD_90(dir_out, h, roughness);
    Real FD = (1. + (FD90 - 1.)*(1. - pow(fabs(dot(norm, dir)),5)));
    return FD;
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

    // Homework 1: implement this!
    //std::cout << "dir_in = " << length(dir_in) << "\n";
    //std::cout << "dir_out = " << length(dir_out) << "\n";

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        
    Vector3 h = (dir_in + dir_out)/length((dir_in + dir_out));
    Real FD_in = FD(frame.n, dir_in, dir_out, h, roughness);
    Real FD_out = FD(frame.n, dir_out, dir_out, h, roughness);
    
    //Spectrum f_base_diffuse = bsdf.base_color * 
    //std::cout << "FD = " << FD << " \n";
    Vector3 baseColor = eval(
            bsdf.base_color, 
            vertex.uv, 
            vertex.uv_screen_size, 
            texture_pool);

    Real fac = (FD_in * FD_out/c_PI) * 
        fabs(dot(h, dir_out));

    //if(fac > 1.0)
    //    std::cout << "fac = " << fac << "\n";

    return fac * baseColor;
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
