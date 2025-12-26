# ⚛️ 核反应堆毒物仿真 Pro 

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://streamlit.io)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](./LICENSE)

这是一个基于 Python 和 Streamlit 构建的交互式核反应堆物理仿真工具。它主要用于模拟和分析反应堆运行过程中，裂变产物毒物（**氙-135** 和 **钐-149**）的瞬态行为及其对反应性（Reactivity Worth）的影响。

该工具非常适合核工程专业的学生、教师以及对反应堆物理感兴趣的工程师，用于直观地理解“碘坑效应”、“停堆后钐峰”等复杂物理现象。

## ✨ 主要功能 (Features)

* **双核素体系仿真**：
    * **碘-氙体系 (I-135 → Xe-135)**：模拟高吸收截面毒物的动态平衡与著名的“碘坑”效应。
    * **钷-钐体系 (Pm-149 → Sm-149)**：模拟稳定毒物在长期运行及停堆后的积累特性。
* **交互式工况配置**：
    * 提供预设典型工况：**冷态启动 (Startup)**、**满功率停堆 (Shutdown/Iodine Pit)**、**功率台阶变化 (Step Change)**。
    * 支持**自定义运行历史**：用户可自由配置多个时间段的功率水平和持续时间。
* **实时数值解算**：
    * 基于 `scipy.integrate.odeint` 实时求解 Bateman 微分方程组。
    * 支持自定义初始核素浓度（全零、平衡态或手动输入）。
* **可视化数据分析**：
    * **多维度绘图**：功率历史、核素浓度变化（前体核 vs 毒物核）、反应性价值 (pcm) 曲线。
    * **物理现象标注**：自动检测并标注最大毒物峰值（如碘坑深度与出现时间）。
    * **详细数据表**：支持查看和下载具体的仿真数据点。

## 🛠️ 安装与运行 (Installation & Usage)

### 1. 克隆仓库
```bash
git clone [https://github.com/your-username/reactor-poison-simulation.git](https://github.com/your-username/reactor-poison-simulation.git)
cd reactor-poison-simulation
```

### 2. 创建虚拟环境 (推荐)
```bash
# Windows
python -m venv venv
venv\Scripts\activate

# macOS/Linux
python3 -m venv venv
source venv/bin/activate
```

### 3. 安装依赖
确保目录下包含 `requirements.txt` 文件，然后运行：
```bash
pip install -r requirements.txt
```

### 4. 启动应用
```bash
streamlit run app.py
```
启动后，浏览器将自动打开 http://localhost:8501。

---

## 📦 依赖列表 (Requirements)
请确保 `requirements.txt` 包含以下库：
```text
streamlit
numpy
scipy
matplotlib
```

## 🧠 物理模型简介 (Physics Model)
本程序核心基于 **Bateman 方程** 描述裂变产物的生成与衰变链：

### 1. 碘-氙动力学 (Iodine-Xenon Dynamics)
$$
\begin{aligned}
\frac{dI}{dt} &= \gamma_I \Sigma_f \phi - \lambda_I I \\
\frac{dX}{dt} &= \gamma_X \Sigma_f \phi + \lambda_I I - \lambda_X X - \sigma_a^X \phi X
\end{aligned}
$$

* **Xe-135** 具有极大的热中子吸收截面 ($\approx 2.65 \times 10^6$ barns)。
* **停堆后**：由于 $\sigma_a^X \phi X$ (中子吸收项) 消失，而碘-135 继续衰变生成氙，导致氙浓度先上升后下降，形成 **“碘坑” (Iodine Pit)**。

### 2. 钷-钐动力学 (Promethium-Samarium Dynamics)
$$
\begin{aligned}
\frac{dP}{dt} &= \gamma_P \Sigma_f \phi - \lambda_P P \\
\frac{dS}{dt} &= \lambda_P P - \sigma_a^S \phi S
\end{aligned}
$$

* **Sm-149** 是稳定同位素，不会自发衰变。
* **停堆后**：钐浓度会持续上升到一个更高的稳定值（停堆后钐峰），必须在重新启动后通过中子吸收才能降低。|

## 🤝 贡献 (Contributing)
欢迎提交 Pull Request 或 Issue！如果您发现物理参数有误或有更好的算法优化建议，请随时联系。

## 📄 许可证 (License)
本项目采用 MIT 许可证 
详情请参阅 [LICENSE](LICENSE) 文件

---
Created with ❤️ by [Xinrui Wang]