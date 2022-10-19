#version 330 core
// Colormap fragment shader

in float value;

uniform vec3 colorMapColors[3];

out vec4 color;

void main()
{
    color = vec4(colorMapColors[int(min(2, value * 3))], 1.0F);
}
