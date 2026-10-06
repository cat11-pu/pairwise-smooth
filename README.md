# smooth

一个纯标准库的曲线插值与缓动内核：二次、三次贝塞尔曲线按各自次数的伯恩斯坦权重求值，切分与截取都走德卡斯特
劳细分，弧长用折线弦长逼近并可以反查，因此整条曲线能按等弧长走位与重采样。缓动档位是单位三次贝塞尔时间曲线，
x 坐标用二分法反解成缓动值。全部是浮点运算，不读时钟、不碰随机数、不做任何 I/O。

## 目录

- smooth/core.py 内核实现（Bezier、CurveError、EASING_PRESETS、ease）
- tests/test_core.py 行为测试

## 约定

- 参数一律取 [0, 1]，越界即报错；控制点是 (x, y) 浮点对，一张曲线 2 到 4 个控制点。
- 弧长按细分折线的弦长累加，默认 64 段；at_distance 按弧长反查参数与位置，samples 据此等弧长重采样。
- 缓动档位名与四个控制点见 EASING_PRESETS，取值按 CSS 时间曲线的读法。

## 运行测试

在项目根目录执行：

    python3 -m unittest discover -s tests -v

只使用 Python 标准库，不需要安装任何依赖。
