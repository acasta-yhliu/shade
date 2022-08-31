# SPPM Render

报告参见`docs/report.pdf`

## Build

```
make
./Shade {scene-config}
```

如果你选择使用其他的编译器，请在`Makefile`的第一行进行修改。

## Scene Config

```
<<BEGIN>> <<CONFIG>>
    <<SIZE>> {width} {height}
    <<NRND>> {rounds}
    <<SAMP>> {photon / 1000}
    <<NAME>> {filename}
<<END>>

<<BEGIN>> <<MATERIAL>> {count}
    {materials}
<<END>>

<<BEGIN>> <<SCENE>> {count}
    {objects}
<<END>>

<<BEGIN>> <<CAMERA>>
    <<POSI>> {x} {y} {z}
    <<LOOK>> {x} {y} {z}
    <<LKUP>> {x} {y} {z} 
    <<CONF>> {fov} {aperture} {focus}
<<END>>

<<ENDCONFIG>>
```

