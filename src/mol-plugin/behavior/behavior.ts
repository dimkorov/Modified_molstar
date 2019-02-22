/**
 * Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * @author David Sehnal <david.sehnal@gmail.com>
 */

import { PluginStateTransform, PluginStateObject } from '../state/objects';
import { StateTransformer, StateTransform } from 'mol-state';
import { Task } from 'mol-task';
import { PluginContext } from 'mol-plugin/context';
import { PluginCommand } from '../command';
import { Observable } from 'rxjs';
import { ParamDefinition } from 'mol-util/param-definition';
import { shallowEqual } from 'mol-util';

export { PluginBehavior }

interface PluginBehavior<P = unknown> {
    register(ref: StateTransform.Ref): void,
    unregister(): void,

    /** Update params in place. Optionally return a promise if it depends on an async action. */
    update?(params: P): boolean | Promise<boolean>
}

namespace PluginBehavior {
    export class Root extends PluginStateObject.Create({ name: 'Root', typeClass: 'Root' }) { }
    export class Behavior extends PluginStateObject.CreateBehavior<PluginBehavior>({ name: 'Behavior' }) { }

    export interface Ctor<P = undefined> { new(ctx: PluginContext, params: P): PluginBehavior<P> }

    export interface CreateParams<P> {
        name: string,
        ctor: Ctor<P>,
        canAutoUpdate?: StateTransformer.Definition<Root, Behavior, P>['canAutoUpdate'],
        label?: (params: P) => { label: string, description?: string },
        display: {
            name: string,
            group: string,
            description?: string
        },
        params?(a: Root, globalCtx: PluginContext): { [K in keyof P]: ParamDefinition.Any }
    }

    export function create<P>(params: CreateParams<P>) {
        // TODO: cache groups etc
        return PluginStateTransform.CreateBuiltIn<Root, Behavior, P>({
            name: params.name,
            display: params.display,
            from: [Root],
            to: [Behavior],
            params: params.params,
            apply({ params: p }, ctx: PluginContext) {
                const label = params.label ? params.label(p) : { label: params.display.name, description: params.display.description };
                return new Behavior(new params.ctor(ctx, p), label);
            },
            update({ b, newParams }) {
                return Task.create('Update Behavior', async () => {
                    if (!b.data.update) return StateTransformer.UpdateResult.Unchanged;
                    const updated = await b.data.update(newParams);
                    return updated ? StateTransformer.UpdateResult.Updated : StateTransformer.UpdateResult.Unchanged;
                })
            },
            canAutoUpdate: params.canAutoUpdate
        });
    }

    export function simpleCommandHandler<T>(cmd: PluginCommand<T>, action: (data: T, ctx: PluginContext) => void | Promise<void>) {
        return class implements PluginBehavior<{}> {
            // TODO can't be private due to bug with generating declerations, see https://github.com/Microsoft/TypeScript/issues/17293
            /** private */ sub: PluginCommand.Subscription | undefined = void 0;
            register(): void {
                this.sub = cmd.subscribe(this.ctx, data => action(data, this.ctx));
            }
            unregister(): void {
                if (this.sub) this.sub.unsubscribe();
                this.sub = void 0;
            }
            // TODO can't be private due to bug with generating declerations, see https://github.com/Microsoft/TypeScript/issues/17293
            constructor(/** private */ public ctx: PluginContext) { }
        }
    }

    export abstract class Handler<P = { }> implements PluginBehavior<P> {
        private subs: PluginCommand.Subscription[] = [];
        protected subscribeCommand<T>(cmd: PluginCommand<T>, action: PluginCommand.Action<T>) {
            this.subs.push(cmd.subscribe(this.ctx, action));
        }
        protected subscribeObservable<T>(o: Observable<T>, action: (v: T) => void) {
            this.subs.push(o.subscribe(action));
        }
        protected track<T>(sub: PluginCommand.Subscription) {
            this.subs.push(sub);
        }
        abstract register(): void;
        unregister() {
            for (const s of this.subs) s.unsubscribe();
            this.subs = [];
        }
        update(params: P): boolean | Promise<boolean> {
            if (shallowEqual(params, this.params)) return false;
            this.params = params;
            return true;
        }
        constructor(protected ctx: PluginContext, protected params: P) {
        }
    }
}