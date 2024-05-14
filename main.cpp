#include <Novice.h>

#include <Vector2.h>
#include <Vector3.h>
#include <imgui.h>

#include "MyMath.h"

const char kWindowTitle[] = "LC1B_13_コムロ_リュウヘイ";

struct Mouse {
	Vector2d pos;
	int radius = 8;
	unsigned int color = 0xFFFFFFFF;
};

struct Enemy {
	Vector2 pos;
	Vector2 vel;
	int radius;
	unsigned int color;
	bool isActive;
};

bool IsCollision(const Enemy& enemy, const Mouse& mouse) {

	// 2つの球の中心点間の距離を求める
	float distance = Length(enemy.pos, mouse.pos);
	// 半径の合計よりも短ければ衝突
	if (distance <= enemy.radius + mouse.radius) {
		return true;
	}
	return false;
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = { 0 };
	char preKeys[256] = { 0 };

	Mouse mouse;
	Enemy enemy[5];
	for (int i = 0; i < 5; i++) {
		enemy[i].pos.x = 100.0f;
		enemy[i].pos.y = 100.0f + rand() % 521;
		enemy[i].vel.x = 50.0f + static_cast<float>(rand() % 51);
		enemy[i].radius = 32;
		enemy[i].color = 0xFFFFFFFF;
		enemy[i].isActive = true;
	}

	float currentSpeed = 1.0f;
	float targetSpeed = 1.0f;
	float transitionSpeed = 0.02f; // 速度の変化の速さ

	int score = 0;

	bool isBulletTime = false;
	int bulletTimeTimer = 300;
	int bulletTimeCoolTime = 180;

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		Novice::GetMousePosition(&mouse.pos.x, &mouse.pos.y);

		for (int i = 0; i < 5; i++) {

			// 当たり判定の処理
			if (IsCollision(enemy[i], mouse)) {
				mouse.color = 0xFF0000FF;
				if (Novice::IsTriggerMouse(0)) {
					// 白い敵を打つとオレンジになる
					if (enemy[i].color == 0xFFFFFFFF) {
						enemy[i].color = 0xFFA500FF;
						// オレンジの敵を打つと赤になる
					} else if (enemy[i].color == 0xFFA500FF) {
						enemy[i].color = 0xFF0000FF;
						// 赤い敵を打つと敵のフラグがfalseに
					} else if (enemy[i].color == 0xFF0000FF) {
						enemy[i].isActive = false;
						score = score + 100;
					}
				}
			} else {
				mouse.color = 0xFFFFFFFF;
			}
		}

		// 左矢印キーでバレットタイムを有効
		if (!isBulletTime) {
			if (bulletTimeTimer == 300) {
				if (keys[DIK_LEFT]) {
					isBulletTime = true;
					targetSpeed = 0.1f;
				}
			}
		} else {
			bulletTimeTimer--;
			if (bulletTimeTimer <= 0) {
				isBulletTime = false;
				targetSpeed = 1.0f;
				bulletTimeCoolTime = 180;
			}
		}

		// クールタイムを減らす処理
		if (!isBulletTime && bulletTimeTimer == 0) {
			bulletTimeCoolTime--;
		}

		// バレットタイムとクールタイムの時間を戻す処理
		if (bulletTimeCoolTime <= 0) {

			bulletTimeTimer = 300;
			bulletTimeCoolTime = 180;
		}

		// 徐々にスピードを遅くする・戻す処理
		if (currentSpeed < targetSpeed) {
			currentSpeed += transitionSpeed;
			if (currentSpeed > targetSpeed) {
				currentSpeed = targetSpeed;
			}
		} else if (currentSpeed > targetSpeed) {
			currentSpeed -= transitionSpeed;
			if (currentSpeed < targetSpeed) {
				currentSpeed = targetSpeed;
			}
		}

		for (int i = 0; i < 5; i++) {

			enemy[i].pos.x += enemy[i].vel.x * currentSpeed;

			// 敵が画面外に出ると各値を初期化
			if (enemy[i].pos.x >= 1300.0f) {

				enemy[i].pos.x = -100.0f; // 元の場所に戻す
				enemy[i].pos.y = 100.0f + rand() % 521; // 再びPosをランダムに初期化
				enemy[i].vel.x = 50.0f + static_cast<float>(rand() % 51); // 再びVelをランダムに初期化
				enemy[i].color = 0xFFFFFFFF; // 相手のHP(色)を全快
			}

			// 敵が死ぬと各値を初期化
			if (!enemy[i].isActive) {

				enemy[i].pos.x = -100.0f;
				enemy[i].pos.y = 100.0f + rand() % 521;
				enemy[i].vel.x = 50.0f + static_cast<float>(rand() % 51);
				enemy[i].color = 0xFFFFFFFF;
				enemy[i].isActive = true;
			}
		}

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		Novice::DrawEllipse(mouse.pos.x, mouse.pos.y, mouse.radius, mouse.radius, 0.0f, mouse.color, kFillModeWireFrame);

		for (int i = 0; i < 5; i++) {

			if (enemy[i].isActive) {
				Novice::DrawEllipse(int(enemy[i].pos.x), int(enemy[i].pos.y), enemy[i].radius, enemy[i].radius, 0.0f, enemy[i].color, kFillModeWireFrame);
			}
		}

		Novice::ScreenPrintf(10, 10, "currentSpeed %3.2f", currentSpeed);

		Novice::ScreenPrintf(10, 30, "mousePosX %d", mouse.pos.x);
		Novice::ScreenPrintf(10, 50, "mousePosY %d", mouse.pos.y);

		Novice::ScreenPrintf(10, 100, "isBulletTime %d", isBulletTime);
		Novice::ScreenPrintf(10, 120, "bulletTimeTimer %d", bulletTimeTimer);
		Novice::ScreenPrintf(10, 140, "bulletTimeCoolTime %d", bulletTimeCoolTime);

		Novice::ScreenPrintf(600, 30, "%d", score);

		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}