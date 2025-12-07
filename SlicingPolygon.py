import pygame
import math
from typing import List, Tuple, Dict, Any

Point = Tuple[float, float]


# ===================== 기본 수학 유틸 =====================

def cross(ax, ay, bx, by) -> float:
    return ax * by - ay * bx


def dot(ax, ay, bx, by) -> float:
    return ax * bx + ay * by


def polygon_area(poly: List[Point]) -> float:
    """Shoelace formula – 너무 작은 조각 필터링용."""
    n = len(poly)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        s += x1 * y2 - x2 * y1
    return abs(s) * 0.5


def polygon_centroid(poly: List[Point]) -> Point:
    """간단한 중심점 (정점 평균)."""
    if not poly:
        return (0.0, 0.0)
    sx = sy = 0.0
    for x, y in poly:
        sx += x
        sy += y
    n = len(poly)
    return (sx / n, sy / n)


def slice_polygon(
    poly: List[Point],
    blade_start: Point,
    blade_end: Point,
    eps: float = 1e-6
) -> List[List[Point]]:

    n = len(poly)
    if n < 3:
        return [poly]

    x0, y0 = blade_start
    x1, y1 = blade_end

    Dx = x1 - x0 
    Dy = y1 - y0
    D_len2 = Dx * Dx + Dy * Dy
    if D_len2 < eps * eps:
        # 칼 움직임이 너무 짧으면 무시
        return [poly]

    # blade_start 지점에서 폴리곤의 각 정점까지의 이동한 x,y 를 구함
    # sy, sx의 벡터를 blade의 방향벡터와 외적하면 blade의 궤적 기준 위치를 판단할 수 있음
    signed = []
    for (vx, vy) in poly:
        sx = vx - x0
        sy = vy - y0
        s = cross(Dx, Dy, sx, sy)
        signed.append(s)

    polyA: List[Point] = []
    polyB: List[Point] = []
    current = polyA
    intersections = 0

    for i in range(n):    #매개변수로 받은 폴리곤의 한 정점과 그 다음 정점을 순회하며 선택함
        j = (i + 1) % n
        vx_i, vy_i = poly[i]
        vx_j, vy_j = poly[j]
        si = signed[i]
        sj = signed[j]    #각 정점을 이은 엣지들의 blade 기준 위치를 판단함

        current.append((vx_i, vy_i))

        if abs(si) < eps:
            si = 0.0
        if abs(sj) < eps:
            sj = 0.0

        if si == 0.0 and sj == 0.0:
            continue

        # 엣지들의 부호를 기준으로 판단하므로 그 곱이 음수라면 blade의 위, 아래에 위치
        if si * sj < 0.0:
            Ex = vx_j - vx_i
            Ey = vy_j - vy_i

            denom = cross(Dx, Dy, Ex, Ey)
            if abs(denom) < eps:
                # 평행/거의 평행
                continue

            # 엣지 내부에서 교차하는 지 확인 
            t = -cross(Dx, Dy, vx_i - x0, vy_i - y0) / denom
            if t < -eps or t > 1.0 + eps:
                continue
            
            # 엣지 내부에서 blade 벡터와 교차하는 좌표
            ix = vx_i + t * Ex
            iy = vy_i + t * Ey

            # blade 시작점 부터 내부 교차점까지의 벡터가 blade start, end 길이 내부에 있는지확인(blade)
            wx = ix - x0
            wy = iy - y0
            u = dot(wx, wy, Dx, Dy) / D_len2
            if u < -eps or u > 1.0 + eps:
                continue
            
            # 교차횟수 1 증가
            intersections += 1

            
            def append_if_not_duplicate(poly_list: List[Point], px: float, py: float):
                if not poly_list:
                    poly_list.append((px, py))
                    return
                lx, ly = poly_list[-1]
                if (lx - px) * (lx - px) + (ly - py) * (ly - py) > eps * eps:
                    poly_list.append((px, py))

            #교차점은 분리된 두 폴리곤에 모두 들어가야하므로 양쪽에 추가
            append_if_not_duplicate(polyA, ix, iy)
            append_if_not_duplicate(polyB, ix, iy)

            # 교차점을 넘은 poly의 정점들을 다른 poly에 담기 위해 현재 리스트 변경
            current = polyB if current is polyA else polyA

    if intersections < 2:
        # 교차가 2회 미만이라면 절단되지 않음
        return [poly]

    return [polyA, polyB]


def create_body(
    poly: List[Point],
    vx: float = 0.0,
    vy: float = 0.0,
    ang_vel: float = 0.0,
    last_sliced_id: int = -1
) -> Dict[str, Any]:
    area = polygon_area(poly)

    min_mass = 50.0
    mass = max(area, min_mass)
    return {
        "points": poly,
        "vel": (vx, vy),
        "ang_vel": ang_vel,
        "mass": mass,
        "last_sliced_id": last_sliced_id,
    }


def apply_slice_to_bodies(
    bodies: List[Dict[str, Any]],
    blade_start: Point,
    blade_end: Point,
    blade_speed: float,    
    stroke_id: int,        
    eps: float = 1e-6
) -> List[Dict[str, Any]]:

    x0, y0 = blade_start
    x1, y1 = blade_end

    Dx = x1 - x0
    Dy = y1 - y0
    length = math.hypot(Dx, Dy)
    if length < eps:
        return bodies

    dirx = Dx / length
    diry = Dy / length

    nx = -diry
    ny = dirx

    base_impulse = blade_speed * 8.0  
    max_linear_speed = 1000.0          

    max_ang_vel = 15.0                 

    inherit_linear = 0.3               
    inherit_angular = 0.3              

    new_bodies: List[Dict[str, Any]] = []

    # 전달받은 폴리곤들을 순회하며 동작
    for body in bodies:
        poly = body["points"]
        old_vx, old_vy = body["vel"]
        old_w = body["ang_vel"]
        last_id = body.get("last_sliced_id", -1)

        # 이번 slice id를 받아와서 이번 slice에 잘렸다면 잘리지 않도록 함
        if last_id == stroke_id:
            new_bodies.append(body)
            continue

        #slice_polygon을 호출하여, blade 점들을 기준으로 폴리곤을 slice함        
        pieces = slice_polygon(poly, blade_start, blade_end)

        #잘리지 않아 리스트의 길이가 1인경우 경우 그대로 new body에 저장
        if len(pieces) == 1:
            new_bodies.append(body)
            continue

        # 잘린 경우 처리
        for piece in pieces:
            if len(piece) < 3:
                continue # 잘린 폴리곤의 정점 개수가 3개 미만이면 도형이 아니므로  new bodies에 추가하지 않음
            area = polygon_area(piece) 
            if area <= 1.0:
                continue # 잘린 폴리곤의 면적이 너무 작을 경우 new bodies에 추가하지 않음

            tmp_body = create_body(piece) # 잘린 폴리곤 정점을 가져와 body 생성
            m = tmp_body["mass"] 
            cx, cy = polygon_centroid(piece)

            # 잘린 폴리곤의 중심의 blade 궤적 기준 위, 아래 위치 판단
            s = cross(Dx, Dy, cx - x0, cy - y0)
            side_sign = 1.0 if s >= 0.0 else -1.0

            # 속도를 기본 충격량을 질량으로 나눠 질량이 작을 수록 더 빨리 돌도록 함
            linear_speed = base_impulse / m
            linear_speed = min(linear_speed, max_linear_speed)

            sep_speed = 200.0 / max(math.sqrt(m), 10.0)  # 가벼운 조각일수록 더 벌어짐

            vx = dirx * linear_speed + nx * sep_speed * side_sign
            vy = diry * linear_speed + ny * sep_speed * side_sign

            # 폴리곤이 이미 움직이고 있었다면, 그 속도를 일부 상속받음
            vx += old_vx * inherit_linear
            vy += old_vy * inherit_linear

            # blade 궤적에서부터 중심까지 거리
            signed_distance = s / length  # 궤적 기준 거리 (±)
            mass_scale = max(m * 0.1, 10.0)
            ang_vel = (blade_speed * 0.003) * (signed_distance / mass_scale) # blade의 속도와 거리, 질량을 고려하여 속도결정
            ang_vel += old_w * inherit_angular

            # 각속도 제한
            if ang_vel > max_ang_vel:
                ang_vel = max_ang_vel
            elif ang_vel < -max_ang_vel:
                ang_vel = -max_ang_vel
            
            # 임시 바디에 위 사항을 적용하고 new_bodies에 추가함
            tmp_body["vel"] = (vx, vy)
            tmp_body["ang_vel"] = ang_vel
            tmp_body["last_sliced_id"] = stroke_id

            new_bodies.append(tmp_body)

    return new_bodies

def project_polygon(poly: List[Point], ax: float, ay: float) -> Tuple[float, float]:
    min_p = 1e9
    max_p = -1e9
    for x, y in poly:
        p = x * ax + y * ay
        if p < min_p:
            min_p = p
        if p > max_p:
            max_p = p
    return min_p, max_p


def sat_poly_poly(poly1: List[Point], poly2: List[Point]) -> Tuple[bool, float, float, float]:
    
    #받아온 폴리곤의 정점 개수로 도형인지 판단함
    if len(poly1) < 3 or len(poly2) < 3:
        return (False, 0.0, 0.0, 0.0)

    #SAT 검사를 위한 축을 저장하는 리스트
    axes: List[Tuple[float, float]] = []

    # 폴리곤 하나의 정점을 순회하며 엣지 벡터를 만들고
    # 엣지 벡터에 수직인 노말벡터를 계산하고 정규화 함
    def add_axes_from(poly: List[Point]):
        n = len(poly)
        for i in range(n):
            x1, y1 = poly[i]
            x2, y2 = poly[(i + 1) % n]
            ex = x2 - x1
            ey = y2 - y1
            if abs(ex) < 1e-6 and abs(ey) < 1e-6:
                continue
            # edge에 수직인 축 (노말)
            nx = -ey
            ny = ex
            length = math.hypot(nx, ny)
            if length < 1e-6:
                continue
            nx /= length
            ny /= length
            axes.append((nx, ny))

    # 충돌 검사할 모든 엣지들에 대하여 검사 축 생성 및 추가
    add_axes_from(poly1)
    add_axes_from(poly2)

    c1x, c1y = polygon_centroid(poly1)
    c2x, c2y = polygon_centroid(poly2)
    dx = c2x - c1x
    dy = c2y - c1y

    min_overlap = 1e9
    mtv_axis = (0.0, 0.0)

    for (ax, ay) in axes:
        #axes에 저장된 검사축으로 각 폴리곤들을 투영
        # poly1 투영
        min1, max1 = project_polygon(poly1, ax, ay)
        # poly2 투영
        min2, max2 = project_polygon(poly2, ax, ay)

        # 분리축 체크
        # 1의 최대보다 2의 최소가 작거나 2의 최대보다 1의 최소가 큰 경우 false 반환
        if max1 < min2 or max2 < min1:
            return (False, 0.0, 0.0, 0.0)

        # 겹친 정도 계산
        # 더 작은 max와 더 큰 min 사이의 길이
        overlap = min(max1, max2) - max(min1, min2) 
        if overlap < min_overlap:
            min_overlap = overlap
            # 축 방향을 poly1 -> poly2 방향으로 맞추기
            # 1 -> 2 의 중심으로 의 벡터와 각 축의 내적을 구해서 방향을 1->2가 양수로 통일함
            if dx * ax + dy * ay < 0.0: 
                ax = -ax
                ay = -ay
            mtv_axis = (ax, ay)

    # 충돌 여부와 밀어낼 방향과 깊이 반환
    return (True, mtv_axis[0], mtv_axis[1], min_overlap)


def translate_poly(poly: List[Point], dx: float, dy: float) -> List[Point]:
    # 델타만큼 폴리곤 정점 이동
    return [(x + dx, y + dy) for x, y in poly]


def resolve_collisions_poly(bodies: List[Dict[str, Any]], iterations: int = 2):
    # 조각이 2개 미만이면 return
    if len(bodies) < 2:
        return
    
    for _ in range(iterations):
        for i in range(len(bodies)):
            for j in range(i + 1, len(bodies)):
                # 충돌 검사할 두 조각
                bi = bodies[i]
                bj = bodies[j]

                poly1 = bi["points"]
                poly2 = bj["points"]
                mi = bi["mass"]
                mj = bj["mass"]

                # 함수의 리턴값을 받아 충돌여부, MTV 방향 단위 벡터, 겹친 거리를 저장함
                colliding, nx, ny, depth = sat_poly_poly(poly1, poly2)
                if not colliding or depth <= 0.0:
                    continue

                total_mass = mi + mj
                # 총 질량에 비하 각 도형의 질량을 나눠 작은 질량일 수록 먼 거리를 밀려나도록 설정
                move1 = -depth * (mj / total_mass) * 0.5
                move2 = depth * (mi / total_mass) * 0.5

                #두 조각의 모든 정점을 move 만큼 이동
                bi["points"] = translate_poly(poly1, nx * move1, ny * move1)
                bj["points"] = translate_poly(poly2, nx * move2, ny * move2)

                v1x, v1y = bi["vel"]
                v2x, v2y = bj["vel"]

                # 2기준 상대 속도를 구하여 법선 방향에 투영
                # 음수라면 1->2로 접근중이라는 뜻
                rvx = v2x - v1x
                rvy = v2y - v1y
                rel_vn = rvx * nx + rvy * ny

                if rel_vn < 0.0:  # 서로에게 다가오는 상황에만 충돌시 반사 속도를 추가함
                    restitution = 0.2
                    j_imp = -(1.0 + restitution) * rel_vn
                    j_imp /= (1.0 / mi + 1.0 / mj)

                    imp_x = j_imp * nx
                    imp_y = j_imp * ny

                    v1x -= imp_x / mi
                    v1y -= imp_y / mi
                    v2x += imp_x / mj
                    v2y += imp_y / mj

                    bi["vel"] = (v1x, v1y)
                    bj["vel"] = (v2x, v2y)

def main():
    pygame.init()
    W, H = 800, 600
    screen = pygame.display.set_mode((W, H))
    pygame.display.set_caption("Slicing Polygon, SAT Collision")
    clock = pygame.time.Clock()

    # 초기 하나의 몸체 (볼록 폴리곤)
    initial_poly = [
        (200, 150),
        (500, 120),
        (600, 300),
        (450, 450),
        (220, 380),
    ]
    bodies: List[Dict[str, Any]] = [create_body(initial_poly)]

    slicing = False
    stroke_start: Point = (0, 0)
    prev_mouse_pos: Point = (0, 0)
    curr_mouse_pos: Point = (0, 0)

    trail: List[Point] = []
    max_trail_len = 40

    stroke_id = 0  # 스윙마다 1 증가

    running = True
    while running:
        dt_ms = clock.tick(60)
        dt = dt_ms / 1000.0  # 초 단위

        # 이벤트 처리
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            #마우스 입력 처리 시작 위치 기록
            elif event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 1:
                    slicing = True
                    stroke_id += 1
                    stroke_start = event.pos
                    prev_mouse_pos = event.pos
                    curr_mouse_pos = event.pos
                    trail.clear()
                    trail.append(event.pos)

            # 드래그 중일 때 위치 기록
            elif event.type == pygame.MOUSEMOTION:
                if slicing:
                    curr_mouse_pos = event.pos
                    trail.append(curr_mouse_pos)
                    if len(trail) > max_trail_len:
                        trail.pop(0)

                    # 마우스 속도 계산
                    x0, y0 = prev_mouse_pos
                    x1, y1 = curr_mouse_pos
                    dx = x1 - x0
                    dy = y1 - y0
                    dist = math.hypot(dx, dy)

                    if dist > 1.0 and dt > 0.0:
                        # 프레임 당 이동 거리를 시간으로 나눈 속도
                        blade_speed = dist / dt  # px/s
                        # 마우스 정보를 넣고 궤적에 body와 교차했으면 슬라이스
                        bodies = apply_slice_to_bodies(
                            bodies,
                            stroke_start,
                            curr_mouse_pos,
                            blade_speed,
                            stroke_id
                        )

                    prev_mouse_pos = curr_mouse_pos
            # 마우스 클릭 종료
            elif event.type == pygame.MOUSEBUTTONUP:
                if event.button == 1:
                    slicing = False
                    trail.clear()

        # 속도, 각속도 감쇠 계수
        lin_damping = 0.985
        ang_damping = 0.985
        vel_eps2 = 1e-2

        # 각 바디마다 처리함
        for body in bodies:
            poly = body["points"]
            vx, vy = body["vel"]
            w = body["ang_vel"]

            # 선형 이동
            moved_poly: List[Point] = []
            for (x, y) in poly:
                moved_poly.append((x + vx * dt, y + vy * dt))

            # 회전
            if abs(w) > 1e-4:
                cx, cy = polygon_centroid(moved_poly)
                cos_a = math.cos(w * dt)
                sin_a = math.sin(w * dt)

                rotated_poly: List[Point] = []
                for (x, y) in moved_poly:
                    rx = x - cx
                    ry = y - cy
                    x2 = cx + rx * cos_a - ry * sin_a
                    y2 = cy + rx * sin_a + ry * cos_a
                    rotated_poly.append((x2, y2))
                body["points"] = rotated_poly
            else:
                body["points"] = moved_poly

            # 감쇠
            vx *= lin_damping
            vy *= lin_damping
            w *= ang_damping

            # 너무 작은 속도는 0으로
            if vx * vx + vy * vy < vel_eps2:
                vx = vy = 0.0
            if abs(w) < 1e-3:
                w = 0.0

            body["vel"] = (vx, vy)
            body["ang_vel"] = w

        # 모든 body들에 대해 충돌 검사
        resolve_collisions_poly(bodies, iterations=2)

        # 화면 밖으로 완전히 나간 조각 정리
        def is_off_screen(points: List[Point]) -> bool:
            xs = [p[0] for p in points]
            ys = [p[1] for p in points]
            if max(xs) < -100 or min(xs) > W + 100:
                return True
            if max(ys) < -100 or min(ys) > H + 100:
                return True
            return False

        bodies = [b for b in bodies if not is_off_screen(b["points"])]

        # 그리기
        screen.fill((0, 0, 0))

        colors = [(255, 255, 255)]

        for i, body in enumerate(bodies):
            poly = body["points"]
            if len(poly) >= 3:
                color = colors[i % len(colors)]
                pygame.draw.polygon(screen, color, poly, 0)
                pygame.draw.polygon(screen, (255, 255, 255), poly, 2)

        # 칼 궤적(전체 스윙 경로) 그리기
        if len(trail) > 1:
            pygame.draw.lines(screen, (255, 105, 180), False, trail, 2)

        font = pygame.font.SysFont("consolas", 18)
        text = font.render(
            "drag mouse to slice polygon ",
            True,
            (255, 255, 255),
        )
        screen.blit(text, (10, 10))

        pygame.display.flip()

        keys = pygame.key.get_pressed()
        if keys[pygame.K_ESCAPE]:
            running = False

    pygame.quit()


if __name__ == "__main__":
    main()
