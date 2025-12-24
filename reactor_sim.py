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
# 4. 辅助绘图函数
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
    
    # 自动标注峰值（最负反应性）
    min_idx = np.argmin(worth) 
    peak_val = worth[min_idx]
    peak_time = time_h[min_idx]
    
    # 仅当峰值显著时标注
    if peak_val < -100:
        ax3.annotate(f'Max Poison: {peak_val:.0f} pcm\nT={peak_time:.1f}h', 
                     xy=(peak_time, peak_val), 
                     xytext=(0, 40),            
                     textcoords='offset points', 
                     ha='center',
                     arrowprops=dict(facecolor='black', arrowstyle='->', connectionstyle='arc3')) 

    plt.tight_layout()
    return fig

# ==========================================
# 5. Streamlit 用户界面主程序
# ==========================================
def main():
    st.set_page_config(page_title="核反应堆毒物仿真 Pro", layout="wide", page_icon="⚛️")
    
    st.title("⚛️ 核反应堆裂变产物瞬态分析系统")
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
        
        if init_mode == "平衡态 (基于满功率)":
            phi = FULL_POWER_FLUX
            I_eq = const.GAMMA_I * const.SIGMA_F * phi / const.LAMBDA_I
            X_eq = (const.GAMMA_X + const.GAMMA_I) * const.SIGMA_F * phi / (const.LAMBDA_X + const.SIGMA_A_X * phi)
            P_eq = const.GAMMA_PM * const.SIGMA_F * phi / const.LAMBDA_PM
            S_eq = const.GAMMA_PM * const.SIGMA_F / const.SIGMA_A_S 
            y0 = [I_eq, X_eq, P_eq, S_eq]
            st.success(f"已加载平衡态:\nXe: {X_eq:.2e}\nSm: {S_eq:.2e}")
            
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
    if st.button("🚀 开始计算 (Run Simulation)", type="primary", width="stretch"):
        
        # 1. 运行数值模拟
        t_arr, y_arr = simulate_transient(stages_input, y0)
        
        # 2. 数据处理
        time_hours = t_arr / 3600.0
        I_conc = y_arr[:, 0]
        X_conc = y_arr[:, 1]
        P_conc = y_arr[:, 2]
        S_conc = y_arr[:, 3]
        
        # 3. 反应性价值计算 (估算 pcm)
        Sigma_Xe = X_conc * const.SIGMA_A_X
        Sigma_Sm = S_conc * const.SIGMA_A_S
        
        # 计算满功率参考值用于标定 (假设满功率平衡氙价值约 -2800 pcm)
        phi_ref = FULL_POWER_FLUX
        X_eq_ref = (const.GAMMA_X + const.GAMMA_I) * const.SIGMA_F * phi_ref / (const.LAMBDA_X + const.SIGMA_A_X * phi_ref)
        Sigma_Xe_ref = X_eq_ref * const.SIGMA_A_X
        pcm_scaling = -2800.0 / Sigma_Xe_ref if Sigma_Xe_ref > 0 else 0
        
        Rho_Xe = Sigma_Xe * pcm_scaling
        Rho_Sm = Sigma_Sm * pcm_scaling 

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

            # --- 动态物理解释逻辑 ---
            st.markdown("### 💡 物理现象深度解析")
            if "启动" in scenario:
                st.info("""
                **工况：新堆冷态启动 (Startup)**
                1. **积累过程**：初始时刻 I-135 和 Xe-135 均为 0。随着功率提升，I-135 (前体核) 首先由裂变迅速积累。
                2. **滞后效应**：Xe-135 的积累滞后于 I-135，因为它的主要来源是 I-135 的衰变。
                3. **平衡态**：约 40-50 小时后，Xe-135 的生成（裂变+I衰变）与消失（中子吸收+衰变）达到平衡，反应性价值趋于稳定。
                """)
            elif "停堆" in scenario:
                st.warning(r"""
                **工况：满功率停堆 - 碘坑效应 (Iodine Pit)** 
                1. **消失项归零**：停堆瞬间，中子通量 $\phi \to 0$，Xe-135 的主要消失途径（中子吸收 $\sigma_a \phi X$）立刻停止。
                2. **生成项持续**：堆内积累的大量 I-135 继续以 6.6 小时的半衰期衰变为 Xe-135。
                3. **结果**：生成速率 > 消失速率（仅剩衰变），导致 Xe-135 浓度不降反升，在停堆后 **9-12小时** 出现峰值（即“碘坑”），随后才随时间衰减。
                """)
            elif "台阶" in scenario:
                st.info("""
                **工况：功率台阶变化 (Step Change)**
                1. **瞬态超调**：功率下降瞬间，中子吸收能力减弱，但 I-135 的积累量仍处于高功率水平，导致 Xe-135 浓度短时间内先上升。
                2. **趋向新平衡**：随着 I-135 浓度随裂变率降低而下降，Xe-135 最终会稳定在对应低功率的新平衡点。
                """)
            else:
                st.info("当前为自定义输入模式，请观察曲线中生成项（I衰变）与消失项（吸收+衰变）的竞争关系。")

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

            # --- 动态物理解释逻辑 ---
            st.markdown("### 💡 物理现象深度解析")
            if "启动" in scenario:
                st.info("""
                **工况：新堆冷态启动**
                1. **长周期积累**：Pm-149 和 Sm-149 达到平衡的时间比碘-氙体系长得多（数周时间）。
                2. **平衡特性**：Sm-149 是稳定核素（不衰变）。值得注意的是，**平衡钐浓度与中子通量大小无关**，仅取决于核截面参数。
                """)
            elif "停堆" in scenario:
                st.error(r"""
                **工况：满功率停堆 - 停堆后钐峰**
                1. **只增不减**：停堆后，Sm-149 的中子吸收项消失（$\sigma_a \phi S = 0$），即不再被“烧掉”。
                2. **持续积累**：Pm-149 继续衰变补充 Sm-149，导致 Sm-149 浓度上升至比运行水平更高的峰值。
                3. **永久性**：与 Xe-135 不同，Sm-149 是稳定的。除非重新启动反应堆将其烧掉，否则**它将永久维持在这个高浓度水平**。
                """)
            elif "台阶" in scenario:
                st.info("""
                **工况：功率台阶变化**
                1. **过渡过程**：功率下降初期，由于 Pm-149 的积累存量，Sm-149 浓度会暂时上升。
                2. **最终状态**：长期来看，Sm-149 会回归到与之前几乎相同的平衡值（因为平衡钐浓度与功率水平无关）。
                """)
            else:
                st.info("自定义模式分析：注意观察 Sm-149 作为稳定毒物的积累特性（只通过中子吸收消失）。")

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