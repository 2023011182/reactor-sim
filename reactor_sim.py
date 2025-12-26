import streamlit as st
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# ==========================================
# 1. 物理常数定义
# ==========================================
class ReactorConstants:
    # 衰变常数 (s^-1)
    LAMBDA_I = 2.9306e-5   # I-135 -> Xe-135
    LAMBDA_X = 2.1065e-5   # Xe-135 -> Cs-135
    LAMBDA_PM = 3.6274e-6  # Pm-149 -> Sm-149
    
    # 裂变产额 (直接产额 + 累积产额修正)
    GAMMA_I = 0.0639    
    GAMMA_X = 0.00237   
    GAMMA_PM = 0.01071  
    
    # 微观吸收截面 (cm^2) -> 1 barn = 1e-24 cm^2
    SIGMA_A_X = 2.65e6 * 1e-24  
    SIGMA_A_S = 4.01e4 * 1e-24  
    
    # 宏观裂变截面 (cm^-1, 假设值)
    SIGMA_F = 0.098 
    
    # [新增] 每次裂变产生的平均中子数 (U-235 热裂变典型值)
    NU = 2.43

# ==========================================
# 2. 物理模型核心：Bateman 方程组
# ==========================================
def poison_derivatives(y, t, phi, const):
    I, X, P, S = y
    
    # --- 碘-氙 体系 (I-135 -> Xe-135) ---
    # dI/dt = 裂变产出 - 衰变
    dIdt = const.GAMMA_I * const.SIGMA_F * phi - const.LAMBDA_I * I
    
    # dX/dt = 裂变产出 + 碘衰变补给 - 自发衰变 - 中子吸收
    dXdt = (const.GAMMA_X * const.SIGMA_F * phi + 
            const.LAMBDA_I * I - 
            const.LAMBDA_X * X - 
            const.SIGMA_A_X * phi * X)
            
    # --- 钷-钐 体系 (Pm-149 -> Sm-149) ---
    # dP/dt = 裂变产出 - 衰变
    dPdt = const.GAMMA_PM * const.SIGMA_F * phi - const.LAMBDA_PM * P
    
    # dS/dt = 钷衰变补给 - 中子吸收 (钐是稳定核素，无自发衰变项)
    dSdt = const.LAMBDA_PM * P - const.SIGMA_A_S * phi * S
    
    return [dIdt, dXdt, dPdt, dSdt]

# ==========================================
# 3. 仿真控制逻辑
# ==========================================
def simulate_transient(power_history, initial_state):
    results = []  
    time_points = []  
    current_state = initial_state
    current_time = 0
    
    const = ReactorConstants()
    
    for duration, flux in power_history:
        if duration <= 0: continue
        
        # 自动调整步长，保证绘图平滑
        steps = int(duration * 20) + 10
        t_span = np.linspace(0, duration * 3600, steps) 
        
        # 求解 ODE
        sol = odeint(poison_derivatives, current_state, t_span, args=(flux, const))
        
        abs_time = t_span + current_time
        
        if len(results) == 0:
            results.append(sol)
            time_points.append(abs_time)
        else:
            # 去掉重复的接续点
            results.append(sol[1:])
            time_points.append(abs_time[1:])
            
        current_state = sol[-1]
        current_time += duration * 3600
        
    return np.concatenate(time_points), np.concatenate(results)

# ==========================================
# 4. 辅助绘图函数 (已修改：同时标注极大值和极小值)
# ==========================================
def plot_system_response(time_h, power_x, power_y, 
                          precursor_conc, daughter_conc, worth, 
                          names, colors):
    # 创建 3行1列 的图表布局
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [1, 2, 2]})
    
    # --- 子图 1: 功率历史 ---
    ax1.set_ylabel('Power (%)', fontweight='bold')
    ax1.fill_between(power_x, power_y, color='gray', alpha=0.2)
    ax1.plot(power_x, power_y, color='black', linewidth=1.5)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.set_ylim(0, 120)
    ax1.set_title(f"{names[1]} Transient Dynamics Analysis", fontsize=14)

    # --- 子图 2: 核素浓度 (双Y轴) ---
    ax2_r = ax2.twinx()
    l1 = ax2.plot(time_h, precursor_conc, color=colors[0], linestyle='--', linewidth=2, label=f'{names[0]} (Precursor)')
    ax2.set_ylabel(f'{names[0]} Conc. (atoms/cm³)', color=colors[0], fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=colors[0])
    
    l2 = ax2_r.plot(time_h, daughter_conc, color=colors[1], linewidth=2.5, label=f'{names[1]} (Poison)')
    ax2_r.set_ylabel(f'{names[1]} Conc. (atoms/cm³)', color=colors[1], fontweight='bold')
    ax2_r.tick_params(axis='y', labelcolor=colors[1])
    
    # 合并图例
    lns = l1 + l2
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='best')
    ax2.grid(True, linestyle='--', alpha=0.6)

    # --- 子图 3: 反应性价值 (pcm) ---
    ax3.plot(time_h, worth, color='darkred', linewidth=2)
    ax3.set_ylabel('Reactivity (pcm)', color='darkred', fontweight='bold')
    ax3.set_xlabel('Time (Hours)', fontsize=12)    
    ax3.fill_between(time_h, worth, 0, color='darkred', alpha=0.1)
    ax3.grid(True, linestyle='--', alpha=0.6)
    
    # ================= 修改部分：标注最大值与最小值 =================
    
    # 1. 寻找最小值 (最负反应性，即毒物峰值)
    min_idx = np.argmin(worth) 
    min_val = worth[min_idx]
    min_time = time_h[min_idx]
    
    # 标注最小值 (Text 上浮)
    ax3.annotate(f'Min: {min_val:.0f}\nT={min_time:.1f}h', 
                 xy=(min_time, min_val), 
                 xytext=(0, 40),            
                 textcoords='offset points', 
                 ha='center', va='bottom',
                 arrowprops=dict(facecolor='blue', arrowstyle='->', connectionstyle='arc3'),
                 fontsize=9, color='blue', fontweight='bold') 

    # 2. 寻找最大值 (最正/接近零反应性，即毒物低谷)
    max_idx = np.argmax(worth)
    max_val = worth[max_idx]
    max_time = time_h[max_idx]

    # 标注最大值 (Text 下沉)，仅当最大值点与最小值点不同时标注
    if max_idx != min_idx:
        ax3.annotate(f'Max: {max_val:.0f}\nT={max_time:.1f}h', 
                     xy=(max_time, max_val), 
                     xytext=(0, -40),            
                     textcoords='offset points', 
                     ha='center', va='top',
                     arrowprops=dict(facecolor='red', arrowstyle='->', connectionstyle='arc3'),
                     fontsize=9, color='red', fontweight='bold')
    # =============================================================

    plt.tight_layout()
    return fig

# ==========================================
# 5. Streamlit 用户界面主程序
# ==========================================
def main():
    st.set_page_config(page_title="核反应堆毒物仿真 Pro", layout="wide", page_icon="⚛️")
    
    st.title("⚛️ 核反应堆毒物仿真系统")
    st.markdown("---")
    
    # --- Session State 辅助函数 (用于滑块与数字框同步) ---
    def update_slider(key_slider, key_num):
        st.session_state[key_slider] = st.session_state[key_num]

    def update_num(key_slider, key_num):
        st.session_state[key_num] = st.session_state[key_slider]

    # --- 侧边栏配置区 ---
    with st.sidebar:
        st.header("1. 场景与工况")
        
        scenario = st.selectbox(
            "选择预设典型工况:",
            ["自定义输入", "新堆冷态启动 (Startup)", "满功率停堆-碘坑 (Shutdown)", "功率台阶变化 (Step Change)"]
        )

        st.subheader("2. 堆芯参数")
        FULL_POWER_FLUX = st.number_input(
            "满功率热中子通量", 
            value=3.0e13, 
            format="%.1e",
            step=0.1e13
        )
        
        st.subheader("3. 运行历史配置")
        default_stages = []
        
        # 根据选择的工况预设参数
        if scenario == "新堆冷态启动 (Startup)":
            default_stages = [(100.0, 100)] 
        elif scenario == "满功率停堆-碘坑 (Shutdown)":
            default_stages = [(50.0, 100), (50.0, 0)]
        elif scenario == "功率台阶变化 (Step Change)":
            default_stages = [(40.0, 100), (40.0, 50), (40.0, 100)]
        else: # 自定义默认
            default_stages = [(50.0, 100), (24.0, 0)]

        num_stages = st.number_input("阶段数量", min_value=1, max_value=10, value=len(default_stages))
        stages_input = []
        
        # 动态生成输入框
        for i in range(num_stages):
            st.markdown(f"**阶段 {i+1}**")
            col1, col2 = st.columns([1, 1.2]) 
            
            #以此阶段的预设值为初值，防止越界
            def_dur = default_stages[i][0] if i < len(default_stages) else 10.0
            def_pow = int(default_stages[i][1]) if i < len(default_stages) else 0
            
            key_slider = f"slider_p{i}"
            key_num = f"num_p{i}"
            
            # 初始化 session_state
            if key_slider not in st.session_state:
                st.session_state[key_slider] = def_pow
            if key_num not in st.session_state:
                st.session_state[key_num] = def_pow

            with col1:
                dur = st.number_input(
                    f"时长(h)", 
                    value=float(def_dur), 
                    min_value=0.1, 
                    step=1.0, 
                    key=f"d{i}"
                )
            
            with col2:
                # 数字输入框
                st.number_input(
                    f"功率(%)", 
                    min_value=0, max_value=120, 
                    key=key_num,
                    on_change=update_slider, 
                    args=(key_slider, key_num)
                )
                # 滑块
                st.slider(
                    "功率调节", 
                    min_value=0, max_value=120, 
                    key=key_slider,
                    on_change=update_num,
                    args=(key_slider, key_num),
                    label_visibility="collapsed"
                )
            
            current_p = st.session_state[key_num]
            stages_input.append((dur, FULL_POWER_FLUX * (current_p / 100.0)))
            st.divider()

        st.subheader("4. 初始核素浓度")
        init_mode = st.radio("选择初始条件", ["新堆芯 (全零)", "平衡态 (基于满功率)", "自定义数值"])
        
        y0 = [0.0, 0.0, 0.0, 0.0]
        const = ReactorConstants()
        
        # 预先计算满功率平衡值（用于后续校准和初始值）
        phi_ref = FULL_POWER_FLUX
        I_eq_ref = const.GAMMA_I * const.SIGMA_F * phi_ref / const.LAMBDA_I
        X_eq_ref = (const.GAMMA_X + const.GAMMA_I) * const.SIGMA_F * phi_ref / (const.LAMBDA_X + const.SIGMA_A_X * phi_ref)
        P_eq_ref = const.GAMMA_PM * const.SIGMA_F * phi_ref / const.LAMBDA_PM
        S_eq_ref = const.GAMMA_PM * const.SIGMA_F / const.SIGMA_A_S
        
        if init_mode == "平衡态 (基于满功率)":
            y0 = [I_eq_ref, X_eq_ref, P_eq_ref, S_eq_ref]
            st.success(f"已加载平衡态:\nXe: {X_eq_ref:.2e}\nSm: {S_eq_ref:.2e}")
            
        elif init_mode == "自定义数值":
            st.markdown("初始原子数密度 (atoms/cm³):")
            c1, c2 = st.columns(2)
            y0 = [
                c1.number_input("I-135", 0.0, format="%.1e"),
                c2.number_input("Xe-135", 0.0, format="%.1e"),
                c1.number_input("Pm-149", 0.0, format="%.1e"),
                c2.number_input("Sm-149", 0.0, format="%.1e")
            ]

    # --- 主界面：运行与结果展示 ---
    if st.button("🚀 开始计算", type="primary", width="stretch"):
        
        # 1. 运行数值模拟
        t_arr, y_arr = simulate_transient(stages_input, y0)
        
        # 2. 数据处理
        time_hours = t_arr / 3600.0
        I_conc = y_arr[:, 0]
        X_conc = y_arr[:, 1]
        P_conc = y_arr[:, 2]
        S_conc = y_arr[:, 3]
        
        # 3. 反应性价值计算 
        
        # (a) 计算中子产生项 (Macroscopic Production Cross Section)
        Macroscopic_Production = const.NU * const.SIGMA_F
        
        # (b) 计算基准状态下的本底吸收 (Calibration)
        Sigma_Xe_ref_abs = X_eq_ref * const.SIGMA_A_X
        Sigma_Sm_ref_abs = S_eq_ref * const.SIGMA_A_S
        Sigma_structure = Macroscopic_Production - (Sigma_Xe_ref_abs + Sigma_Sm_ref_abs)
        
        # (c) 计算瞬态过程中的 k_eff 和 反应性
        Sigma_Xe_t = X_conc * const.SIGMA_A_X
        Sigma_Sm_t = S_conc * const.SIGMA_A_S
        
        # 计算针对单一毒物的 k_eff (为了绘图解耦，采用控制变量法)
        k_Xe_only = Macroscopic_Production / (Sigma_structure + Sigma_Xe_t + Sigma_Sm_ref_abs)
        k_Sm_only = Macroscopic_Production / (Sigma_structure + Sigma_Xe_ref_abs + Sigma_Sm_t)
        
        # 转换为反应性 rho = (k-1)/k * 1e5 (pcm)
        Rho_Xe = ((k_Xe_only - 1.0) / k_Xe_only) * 1e5
        Rho_Sm = ((k_Sm_only - 1.0) / k_Sm_only) * 1e5

        # 4. 构建功率绘图数据 (使其为台阶状)
        power_x = [0]
        power_y = [stages_input[0][1]/FULL_POWER_FLUX*100]
        curr_t = 0
        for dur, flx in stages_input:
            p = flx / FULL_POWER_FLUX * 100
            power_x.extend([curr_t, curr_t + dur])
            power_y.extend([p, p])
            curr_t += dur

        # 5. 结果展示 Tab 页
        tab1, tab2, tab3 = st.tabs(["📉 碘-氙 (I-Xe) 动态", "📉 钷-钐 (Pm-Sm) 动态", "📋 详细数据表"])
        
        # ==================================
        # Tab 1: 碘-氙 系统及动态物理解释
        # ==================================
        with tab1:
            st.markdown("#### 碘-135 衰变至 氙-135 过程分析")
            
            fig1 = plot_system_response(
                time_hours, power_x, power_y, 
                I_conc, X_conc, Rho_Xe,
                names=["Iodine-135", "Xenon-135"],
                colors=["tab:orange", "tab:red"]
            )
            st.pyplot(fig1)

            # --- 动态物理解释逻辑---
            st.markdown("### 💡 物理现象深度解析")
            
            if "启动" in scenario:
                st.info(r"""
                **工况：新堆冷态启动 (Startup)**
                1.  **前体核主导 (T < 10h)**：启动初期，I-135 直接由裂变产生 ($\gamma_I \Sigma_f \phi$)，积累迅速。而 Xe-135 主要源于 I-135 的衰变，因此 Xe 的积累表现出明显的 **滞后性 (Lag)**。
                2.  **趋向平衡**：随着 I-135 浓度饱和，其衰变补给率达到最大。Xe-135 的生成项与其消失项（主要是中子吸收 $\sigma_a \phi X$）逐渐抗衡。
                3.  **饱和效应**：最终达到 $dI/dt=0$ 和 $dX/dt=0$ 的动态平衡点，反应性固定在平衡氙毒水平。
                """)
            
            elif "停堆" in scenario:
                st.warning(r"""
                **工况：满功率停堆 - 碘坑效应 (Iodine Pit)**
                1.  **消失机制崩溃**：停堆瞬间，通量 $\phi \to 0$，Xe-135 的主导消失项（中子吸收，占运行时 >90%）立刻归零。仅剩缓慢的自发衰变（$T_{1/2}=9.1h$）。
                2.  **生成项持续**：堆内积累的大量 I-135 继续以 $T_{1/2}=6.6h$ 衰变为 Xe-135。
                3.  **净增量 (Min点)**：由于 **来源(I衰变) > 去路(Xe衰变)**，Xe-135 浓度不降反升，在停堆后 **9-11小时** 达到峰值（死时间窗口），随后因 I-135 耗尽才开始下降。
                """)
            
            elif "台阶" in scenario:
                st.info(r"""
                **工况：功率台阶变化 (Step Change)**
                1.  **降功率**：通量减少导致 Xe 的“燃烧”能力减少，而 I-135 的存量仍处于高位。生成 > 消耗，导致 Xe 浓度暂时上升（**负反应性引入**）。
                2.  **升功率**：通量增加导致 Xe 被剧烈“燃烧”。消耗 > 生成，导致 Xe 浓度急剧下降（**正反应性引入**），需防范瞬态超调。
                """)
            
            else:
                st.info("自定义模式分析：请重点观察曲线中【裂变产出/衰变补给】与【中子吸收/自发衰变】的消长关系。")

        # ==================================
        # Tab 2: 钷-钐 系统及动态物理解释
        # ==================================
        with tab2:
            st.markdown("#### 钷-149 衰变至 钐-149 过程分析")
            
            fig2 = plot_system_response(
                time_hours, power_x, power_y, 
                P_conc, S_conc, Rho_Sm,
                names=["Promethium-149", "Samarium-149"],
                colors=["tab:purple", "tab:blue"]
            )
            st.pyplot(fig2)

            # --- 动态物理解释逻辑 (已修复：使用 r"" 原生字符串) ---
            st.markdown("### 💡 物理现象深度解析")
            
            if "启动" in scenario:
                st.info(r"""
                **工况：新堆冷态启动**
                1.  **长周期积累**：由于 Pm-149 半衰期较长 (53h)，Sm-149 达到平衡需约 3 周时间。
                2.  **平衡特性**：Sm-149 是稳定核素。根据平衡方程，其平衡浓度 $S_{eq}$ 仅取决于裂变产额和截面比值，**与中子通量大小无关**。
                """)
            
            elif "停堆" in scenario:
                st.error(r"""
                **工况：满功率停堆 - 永久中毒**
                1.  **只增不减**：停堆后，Sm-149 的唯一消失途径（中子吸收 $\sigma_a \phi S$）消失。
                2.  **峰值机理**：现存的 Pm-149 全部衰变为 Sm-149，导致 Sm 浓度上升至比运行水平更高的峰值。
                3.  **永久性**：由于 Sm-149 不衰变，该高浓度将**永久保持**，直到下次开堆产生中子将其“烧掉”。
                """)
            
            elif "台阶" in scenario:
                st.info(r"""
                **工况：功率台阶变化 (与 Xe 相反)**
                1.  **降功率**：中子吸收减弱，Pm 衰变补给占优，导致 Sm 浓度**上升**。
                2.  **升功率**：中子吸收增强，加速消耗 Sm，导致 Sm 浓度**下降**。
                *注：Sm 的瞬态变化极其缓慢，通常被视为长期反应性亏损的变化。*
                """)
            
            else:
                st.info("自定义模式分析：注意观察 Sm-149 作为【稳定毒物】的积分积累特性。")
        # ==================================
        # Tab 3: 原始数据
        # ==================================
        with tab3:
            st.dataframe({
                "Time (h)": time_hours,
                "Power (%)": np.interp(time_hours, power_x, power_y),
                "I-135 (at/cm3)": I_conc, 
                "Xe-135 (at/cm3)": X_conc, 
                "Xe Worth (pcm)": Rho_Xe,
                "Pm-149 (at/cm3)": P_conc, 
                "Sm-149 (at/cm3)": S_conc, 
                "Sm Worth (pcm)": Rho_Sm
            }, width="stretch")

if __name__ == "__main__":
    main()