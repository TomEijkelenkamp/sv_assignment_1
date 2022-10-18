#version 330 core
// Height plot fragment shader

in float value;
in float shading;
in float heightChange;

uniform sampler1D textureSampler;

out vec4 color;

void main()
{
    // Grey part scaled by inverse magnitude
    color  = vec4(0.3F, 0.3F, 0.3F, 1.0F) * (1 - heightChange);
    // Color part scaled by magnitude
    color += texture(textureSampler, value) * heightChange;
    // Add shading
    color *= shading;
}
